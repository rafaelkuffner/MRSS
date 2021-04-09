
/*
* Copyright (c) 2021 Victoria University of Wellington CMIC
* @author Ben Allen
* 
*/

#pragma once

#ifndef IOBUFFER_HPP
#define IOBUFFER_HPP

#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <array>
#include <vector>
#include <tuple>
#include <algorithm>
#include <string_view>
#include <charconv>
#include <utility>
#include <limits>
#include <type_traits>
#include <filesystem>

#ifdef _MSC_VER
#define IOB_NOINLINE __declspec(noinline)
#define IOB_FORCEINLINE __forceinline
#else
#define IOB_NOINLINE
#define IOB_FORCEINLINE
#endif

namespace iob {

	enum class endian : unsigned char {
		little, big
	};

	endian native_integer_endian();
	endian native_float_endian();

	namespace detail {

		template <size_t N>
		struct byteswap_impl {
			using type = std::array<unsigned char, N>;
			static void apply(type &x) {
				std::reverse(x.begin(), x.end());
			}
		};

		template <>
		struct byteswap_impl<1> {
			using type = unsigned char;
			static void apply(type &) {}
		};

#if defined(_MSC_VER)
		template <>
		struct byteswap_impl<sizeof(unsigned short)> {
			using type = unsigned short;
			static void apply(type &x) {
				x = _byteswap_ushort(x);
			}
		};

		template <>
		struct byteswap_impl<sizeof(unsigned long)> {
			using type = unsigned long;
			static void apply(type &x) {
				x = _byteswap_ulong(x);
			}
		};

		template <>
		struct byteswap_impl<sizeof(unsigned long long)> {
			using type = unsigned long long;
			static void apply(type &x) {
				x = _byteswap_uint64(x);
			}
		};
#elif defined(__GNUC__)
		template <>
		struct byteswap_impl<sizeof(uint16_t)> {
			using type = uint16_t;
			static void apply(type &x) {
				x = __builtin_bswap16(x);
			}
		};

		template <>
		struct byteswap_impl<sizeof(uint32_t)> {
			using type = uint32_t;
			static void apply(type &x) {
				x = __builtin_bswap32(x);
			}
		};

		template <>
		struct byteswap_impl<sizeof(uint64_t)> {
			using type = uint64_t;
			static void apply(type &x) {
				x = __builtin_bswap64(x);
			}
		};
#endif

		template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
		inline T byteswap(T val) {
			if constexpr (sizeof(T) <= 1) return val;
			using impl = byteswap_impl<sizeof(T)>;
			using swap_t = typename impl::type;
			swap_t x;
			std::memcpy(&x, &val, sizeof(val));
			impl::apply(x);
			std::memcpy(&val, &x, sizeof(val));
			return val;
		}
	}

	class iobuffer {
		// underlying stream is expected to always be in binary mode.
		// append modes should work but will need testing.
	public:
		using uchar = unsigned char;
		using uchar_view = std::basic_string_view<uchar>;
		static constexpr size_t default_bufsize = 131072;
		enum class seek_origin {
			set = SEEK_SET, cur = SEEK_CUR, end = SEEK_END
		};

	private:
		std::vector<uchar> m_buf;
		std::streamoff m_begin_write = 0, m_cursor = 0, m_end_read = 0;
		bool m_eof = false, m_bad = false;
		bool m_reading = false, m_writing = false;
		bool m_report_eof = false;
#ifndef _NDEBUG
		std::streamsize m_put_len = 0;
#endif

	public:
		virtual ~iobuffer() = default;

		iobuffer() = default;

		iobuffer(const iobuffer &) = delete;
		iobuffer & operator=(const iobuffer &) = delete;

		iobuffer(iobuffer &&) = default;
		iobuffer & operator=(iobuffer &&) = default;

		explicit operator bool() const noexcept {
			return !bad() && m_buf.size();
		}

		bool good() const noexcept {
			return !eof() && !bad();
		}

		bool bad() const noexcept {
			return m_bad;
		}

		bool eof() const noexcept {
			// should only report eof after a get that returned less than requested.
			// ie a get(1) that returns the last byte should not cause eof status.
			return m_report_eof;
		}

		IOB_NOINLINE
		std::streampos tell();

		IOB_NOINLINE
		void seek(std::streamoff i, seek_origin origin);

		std::streamsize gavail() const noexcept;

		std::streamsize pavail() const noexcept;

		uchar_view peek(std::streamsize n, std::streamoff base = 0);

		uchar_view get(std::streamsize n);

		std::streamsize get(uchar *buf, std::streamsize n);

		void skip_peeked(uchar_view s);

		void put(uchar_view s);

		uchar * put_begin(std::streamsize n);

		void put_end(std::streamsize n);

		IOB_NOINLINE
		void flush();

	protected:
		void init_open();

	private:
		void read_mode();

		void write_mode();

		IOB_NOINLINE
		void reserve(std::streamsize cap);

		IOB_NOINLINE
		void fill();

		IOB_NOINLINE
		void compact() noexcept;

		std::streamsize capacity() const noexcept;

		virtual bool error_impl() = 0;

		virtual void seek_impl(std::streamoff i, seek_origin origin) = 0;

		virtual std::streampos tell_impl() = 0;

		// if <n bytes are read, eof is assumed
		virtual std::streamsize get_impl(uchar *buf, std::streamsize n) = 0;

		// if <n bytes are written, eof is assumed
		// this should always flush
		virtual std::streamsize put_impl(const uchar *buf, std::streamsize n) = 0;

	};

	// i/o with std::streambuf for compatibility
	class stream_buffer : public iobuffer {
	private:
		// non-owning
		std::streambuf *m_buf = nullptr;
		bool m_bad = false;

	public:
		virtual ~stream_buffer() override;

		stream_buffer() {}

		stream_buffer(stream_buffer &&) noexcept;
		stream_buffer & operator=(stream_buffer &&) noexcept;

		explicit stream_buffer(std::streambuf *);

	private:
		virtual bool error_impl() override;
		virtual void seek_impl(std::streamoff i, seek_origin origin) override;
		virtual std::streampos tell_impl() override;
		virtual std::streamsize get_impl(uchar *buf, std::streamsize n) override;
		virtual std::streamsize put_impl(const uchar *buf, std::streamsize n) override;
	};

	// file i/o with cstdlib
	class file_buffer : public iobuffer {
	public:
#ifdef _WIN32
		using openmode = std::wstring_view;
		static constexpr openmode read = L"rb";
		static constexpr openmode write = L"wb";
		static constexpr openmode append = L"ab";
		static constexpr openmode read_extended = L"rb+";
		static constexpr openmode write_extended = L"wb+";
		static constexpr openmode append_extended = L"ab+";
#else
		using openmode = std::string_view;
		static constexpr openmode read = "rb";
		static constexpr openmode write = "wb";
		static constexpr openmode append = "ab";
		static constexpr openmode read_extended = "rb+";
		static constexpr openmode write_extended = "wb+";
		static constexpr openmode append_extended = "ab+";
#endif
	private:
		FILE *m_fp = nullptr;

	public:
		virtual ~file_buffer() override;

		file_buffer() = default;

		file_buffer(file_buffer &&) noexcept;
		file_buffer & operator=(file_buffer &&) noexcept;

		explicit file_buffer(const std::filesystem::path &fpath, openmode mode);

		bool is_open() const noexcept {
			return bool(m_fp);
		}

	private:
		virtual bool error_impl() override;
		virtual void seek_impl(std::streamoff i, seek_origin origin) override;
		virtual std::streampos tell_impl() override;
		virtual std::streamsize get_impl(uchar *buf, std::streamsize n) override;
		virtual std::streamsize put_impl(const uchar *buf, std::streamsize n) override;
	};

	class text_reader {
	private:
		iobuffer *m_iobuf = nullptr;

	public:
		static constexpr std::string_view whitespace = " \t\r\n";

		text_reader() = default;

		explicit text_reader(iobuffer *iobuf_) : m_iobuf(iobuf_) {}

		iobuffer * iobuf() const noexcept {
			return m_iobuf;
		}

		explicit operator bool() const noexcept {
			return m_iobuf && *m_iobuf;
		}

		bool good() const noexcept {
			return m_iobuf && m_iobuf->good();
		}

		bool bad() const noexcept {
			return !m_iobuf || m_iobuf->bad();
		}

		bool eof() const noexcept {
			return m_iobuf && m_iobuf->eof();
		}
		
		void seek(std::streamoff i, iobuffer::seek_origin origin);

		std::streampos tell() const;

		template <typename T, typename = std::enable_if_t<std::is_unsigned_v<T>>>
		std::string_view peek(T n, std::streamoff base = 0) {
			return peek(std::streamoff(n), base);
		}

		template <typename T, typename = std::enable_if_t<std::is_unsigned_v<T>>>
		std::string_view get(T n) {
			return get(std::streamoff(n));
		}

		std::string_view peek(std::streamsize n, std::streamoff base = 0);

		std::string_view get(std::streamsize n);

		void skip_peeked(std::string_view s);

		// peek up to position indicated by find function
		template <
			typename FindDelim,
			typename = std::enable_if_t<std::is_invocable_r_v<
				std::string::size_type, FindDelim, std::string_view, std::string::size_type
			>>
		>
		std::string_view peek_until(FindDelim &&find, std::streamoff base = 0, std::string::size_type extra = 0) {
			assert(m_iobuf);
			size_t z = 256;
			std::string_view v;
			do {
				auto i0 = v.size();
				v = peek(z, base);
				z += z >> 1;
				if (v.size() > i0) {
					if (auto i = find(v, i0); i != std::string::npos) {
						// found delim
						return peek(i + extra, base);
					}
				} else {
					// no size increase => eof
					break;
				}
			} while (m_iobuf->good());
			// didn't find delim
			return v;
		}

		// peek up to first of `delimchars` (or eof). any delimiter is not returned.
		IOB_NOINLINE
		std::string_view peek_until_any(std::string_view delimchars, std::streamoff base = 0);

		// peek up to `delim` (or eof). any delimiter is not returned.
		IOB_NOINLINE
		std::string_view peek_until_any(char delim, std::streamoff base = 0);

		// peek up to substring `delim` (or eof). any delimiter is not returned.
		IOB_NOINLINE
		std::string_view peek_until_seq(std::string_view delim, std::streamoff base = 0);

		// peek up to first not of `tokchars` (or eof). any delimiter is not returned.
		IOB_NOINLINE
		std::string_view peek_while_any(std::string_view tokchars, std::streamoff base = 0);

		// peek up to first not of `tokchar` (or eof). any delimiter is not returned.
		IOB_NOINLINE
		std::string_view peek_while_any(char tokchar, std::streamoff base = 0);

		// read up to first of `delimchars` (or eof). any delimiter is neither consumed nor returned.
		std::string_view get_until_any(std::string_view delimchars);

		// read up to `delim` (or eof). any delimiter is neither consumed nor returned.
		std::string_view get_until_any(char delim);

		// read up to substring `delim` (or eof). any delimiter is neither consumed nor returned.
		std::string_view get_until_seq(std::string_view delim);

		// read up to first not of `tokchars` (or eof). any delimiter is neither consumed nor returned.
		std::string_view get_while_any(std::string_view tokchars);

		// read up to first not of `tokchar` (or eof). any delimiter is neither consumed nor returned.
		std::string_view get_while_any(char tokchar);

		// read up to `delim` (or eof) and then remove any suffix equal to `strip`.
		// at most one delimiter is consumed from the input but not returned.
		std::string_view get_line(char delim = '\n', char strip = '\r');

		// read up to `delim` (or eof) and then consume at most one delimiter.
		void skip_line(char delim = '\n');

		// check next char
		bool has_any(std::string_view tokchars);

		// check next char
		bool has_any(char c);

		// check for char sequence
		bool has_seq(std::string_view tok);

		// maybe get one char
		std::string_view get_any(std::string_view tokchars);

		// maybe get one char
		std::string_view get_any(char c);

		// maybe get char sequence
		std::string_view get_seq(std::string_view tok);

		std::string_view peek_until_ws();

		std::string_view peek_while_ws();

		std::string_view get_until_ws();

		std::string_view get_while_ws();

		void skip_ws();

		template <typename T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
		IOB_NOINLINE
		std::errc get_int(T &val, int base = 10) {
			static_assert(std::is_integral_v<T>, "must get int as int");
			auto s = peek(32);
			auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), val, base);
			get(end - s.data());
			// TODO return io_error on eof?
			return ec;
		}

		template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
		IOB_FORCEINLINE
		std::errc get_int(T &val, int base = 10) {
			// allow get of integer value as float type
			intmax_t x = 0;
			auto ec = get_int(x, base);
			if (ec != std::errc{}) return ec;
			val = T(x);
			return {};
		}

		template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
		IOB_NOINLINE
		std::errc get_float(T &val, std::chars_format fmt = std::chars_format::general) {
			static_assert(std::is_floating_point_v<T>, "must get float as float");
			auto s = peek(32);
			auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), val, fmt);
			get(end - s.data());
			// TODO return io_error on eof?
			return ec;
		}

		template <typename ...Ts, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_ints(const std::tuple<Ts &...> &vals, Delim delimchars = whitespace, int base = 10) {
			return get_ints_impl(vals, std::index_sequence_for<Ts...>(), delimchars, base);
		}

		template <typename ...Ts, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_floats(const std::tuple<Ts &...> &vals, Delim delimchars = whitespace, std::chars_format fmt = std::chars_format::general) {
			return get_floats_impl(vals, std::index_sequence_for<Ts...>(), delimchars, fmt);
		}

	private:
		template <typename ...Ts, size_t ...Is, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_ints_impl(const std::tuple<Ts &...> &vals, std::index_sequence<Is...>, Delim delimchars, int base = 10) {
			std::errc ec{};
			int r = 0;
			(... && (get_while_any(delimchars), (ec = get_int(std::get<Is>(vals), base)), r += int(ec == std::errc{}), ec == std::errc{}));
			return {ec, r};
		}

		template <typename ...Ts, size_t ...Is, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_floats_impl(const std::tuple<Ts &...> &vals, std::index_sequence<Is...>, Delim delimchars, std::chars_format fmt) {
			std::errc ec{};
			int r = 0;
			(... && (get_while_any(delimchars), (ec = get_float(std::get<Is>(vals), fmt)), r += int(ec == std::errc{}), ec == std::errc{}));
			return {ec, r};
		}

	};

	class binary_reader {
	private:
		iobuffer *m_iobuf = nullptr;
		iob::endian m_endian = endian::little;
		iob::endian m_native_integer_endian = native_integer_endian();
		iob::endian m_native_float_endian = native_float_endian();

	public:
		using uchar = iobuffer::uchar;
		using uchar_view = iobuffer::uchar_view;

		binary_reader() = default;

		explicit binary_reader(iobuffer *iobuf_) : m_iobuf(iobuf_) {}

		iobuffer * iobuf() const noexcept {
			return m_iobuf;
		}

		explicit operator bool() const noexcept {
			return m_iobuf && *m_iobuf;
		}

		bool good() const noexcept {
			return m_iobuf && m_iobuf->good();
		}

		bool bad() const noexcept {
			return !m_iobuf || m_iobuf->bad();
		}

		bool eof() const noexcept {
			return m_iobuf && m_iobuf->eof();
		}

		iob::endian endian() const noexcept {
			return m_endian;
		}

		void endian(iob::endian e) noexcept {
			m_endian = e;
		}

		void seek(std::streamoff i, iobuffer::seek_origin origin);

		std::streampos tell() const;

		template <typename T, typename = std::enable_if_t<std::is_unsigned_v<T>>>
		uchar_view peek(T n, std::streamoff base = 0) {
			return peek(std::streamoff(n), base);
		}

		template <typename T, typename = std::enable_if_t<std::is_unsigned_v<T>>>
		uchar_view get(T n) {
			return get(std::streamoff(n));
		}

		uchar_view peek(std::streamsize n, std::streamoff base = 0);

		uchar_view get(std::streamsize n);

		void skip_peeked(uchar_view s);

		std::string_view get_str(std::streamsize n);

		template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
		T get_val() {
			T val{0};
			const auto nr = m_iobuf->get(reinterpret_cast<uchar *>(&val), sizeof(val));
			if constexpr (std::is_integral_v<T>) {
				if (m_native_integer_endian != m_endian) val = detail::byteswap(val);
			} else if constexpr (std::is_floating_point_v<T>) {
				if (m_native_float_endian != m_endian) val = detail::byteswap(val);
			} else {
				static_assert(false, "bad type for get_val");
			}
			return nr < sizeof(val) ? T{0} : val;
		}

		template <typename ...Ts, typename = std::enable_if_t<(... && std::is_arithmetic_v<Ts>)>>
		IOB_FORCEINLINE
		bool get_vals(Ts &...vals) {
			(void)(..., (vals = get_val<Ts>()));
			return !eof();
		}

		template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
		std::errc get_var_uint(T &val, size_t n) {
			using get_t = uintmax_t;
			assert(n <= sizeof(get_t));
			get_t x{0};
			const auto nr = m_iobuf->get(reinterpret_cast<uchar *>(&x), std::streamsize(n));
			if (nr < std::streamsize(n)) return std::errc::io_error;
			if (m_native_integer_endian != m_endian) x = detail::byteswap(x);
			const size_t shift = CHAR_BIT * (sizeof(get_t) - n);
			if (m_endian == endian::big) x >>= shift;
			if (std::is_integral_v<T> && x > get_t(std::numeric_limits<T>::max())) return std::errc::result_out_of_range;
			val = T(x);
			return std::errc{};
		}

		template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
		std::errc get_var_sint(T &val, size_t n) {
			using get_t = intmax_t;
			assert(n <= sizeof(get_t));
			get_t x{0};
			const auto nr = m_iobuf->get(reinterpret_cast<uchar *>(&x), std::streamsize(n));
			if (nr < std::streamsize(n)) return std::errc::io_error;
			if (m_native_integer_endian != m_endian) x = detail::byteswap(x);
			const size_t shift = CHAR_BIT * (sizeof(get_t) - n);
			if (m_endian == endian::little) x <<= shift;
			static_assert((-1 >> 1) == -1, "need arithmetic right shift");
			x >>= shift;
			if (std::is_integral_v<T> && x < get_t(std::numeric_limits<T>::lowest())) return std::errc::result_out_of_range;
			if (std::is_integral_v<T> && x > get_t(std::numeric_limits<T>::max())) return std::errc::result_out_of_range;
			val = T(x);
			return std::errc{};
		}

		template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
		std::errc get_var_float(T &val, size_t n) {
			static_assert(std::is_floating_point_v<T>, "must get float as float");
			if (n == sizeof(float)) {
				using get_t = float;
				auto x = get_val<get_t>();
				if (eof()) return std::errc::io_error;
				val = T(x);
				return std::errc{};
			} else if (n == sizeof(double)) {
				using get_t = double;
				auto x = get_val<get_t>();
				if (eof()) return std::errc::io_error;
				val = T(x);
				return std::errc{};
			} else {
				assert(false && "bad float size");
				return std::errc::invalid_argument;
			}
		}

	};

	class text_writer {
	private:
		iobuffer *m_iobuf = nullptr;

	public:
		text_writer() = default;

		explicit text_writer(iobuffer *iobuf_) : m_iobuf(iobuf_) {}

		iobuffer * iobuf() const noexcept {
			return m_iobuf;
		}

		explicit operator bool() const noexcept {
			return m_iobuf && *m_iobuf;
		}

		bool good() const noexcept {
			return m_iobuf && m_iobuf->good();
		}

		bool bad() const noexcept {
			return !m_iobuf || m_iobuf->bad();
		}

		void seek(std::streamoff i, iobuffer::seek_origin origin);

		std::streampos tell() const;

		void flush();

		void put(char c);

		void put(std::string_view s);

		char * put_begin(std::streamsize n);

		void put_end(std::streamsize n);

		template <typename T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
		IOB_NOINLINE
		std::errc put_int(const T &val, int base = 10) {
			static_assert(std::is_integral_v<T>, "must put int from int");
			// TODO return io_error if bad?
			char *begin = put_begin(32);
			auto [end, ec] = std::to_chars(begin, begin + 32, val, base);
			// note: on error, to_chars returns ptr == last
			std::streamsize n = (ec == std::errc{}) ? (end - begin) : 0;
			put_end(n);
			return ec;
		}

		template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
		IOB_NOINLINE
		std::errc put_float(const T &val, std::chars_format fmt = std::chars_format::general) {
			static_assert(std::is_floating_point_v<T>, "must put float from float");
			// TODO return io_error if bad?
			char *begin = put_begin(32);
			auto [end, ec] = std::to_chars(begin, begin + 32, val, fmt);
			// note: on error, to_chars returns ptr == last
			std::streamsize n = (ec == std::errc{}) ? (end - begin) : 0;
			put_end(n);
			return ec;
		}

		template <typename T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
		IOB_FORCEINLINE
		std::errc put_float(const T &val, std::chars_format fmt = std::chars_format::general) {
			// allow put of integer type as float representation
			return put_float(double(val), fmt);
		}

		template <typename ...Ts, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> put_ints(const std::tuple<Ts...> &vals, Delim delimseq = ' ', int base = 10) {
			return put_ints_impl(vals, std::index_sequence_for<Ts...>(), delimseq, base);
		}

		template <typename ...Ts, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> put_floats(const std::tuple<Ts...> &vals, Delim delimseq = ' ', std::chars_format fmt = std::chars_format::general) {
			return put_floats_impl(vals, std::index_sequence_for<Ts...>(), delimseq, fmt);
		}

	private:
		template <typename T0, typename ...Ts, size_t I0, size_t ...Is, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> put_ints_impl(const std::tuple<T0, Ts...> &vals, std::index_sequence<I0, Is...>, Delim delimseq, int base = 10) {
			std::errc ec = put_int(std::get<I0>(vals), base);
			if (ec != std::errc{}) return {ec, 0};
			int r = 1;
			(... && (put(delimseq), (ec = put_int(std::get<Is>(vals), base)), r += int(ec == std::errc{}), ec == std::errc{}));
			return {ec, r};
		}

		template <typename T0, typename ...Ts, size_t I0, size_t ...Is, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> put_floats_impl(const std::tuple<T0, Ts...> &vals, std::index_sequence<I0, Is...>, Delim delimseq, std::chars_format fmt) {
			std::errc ec = put_float(std::get<I0>(vals), fmt);
			if (ec != std::errc{}) return {ec, 0};
			int r = 1;
			(... && (put(delimseq), (ec = put_float(std::get<Is>(vals), fmt)), r += int(ec == std::errc{}), ec == std::errc{}));
			return {ec, r};
		}
	};

	constexpr std::string_view rstrip_any(std::string_view s, std::string_view schars = text_reader::whitespace) {
		auto i = s.find_last_not_of(schars);
		if (i == std::string::npos) return {};
		return s.substr(0, i + 1);
	}

	constexpr std::string_view lstrip_any(std::string_view s, std::string_view schars = text_reader::whitespace) {
		auto i = s.find_first_not_of(schars);
		if (i == std::string::npos) return {};
		return s.substr(i);
	}

	constexpr std::string_view strip_any(std::string_view s, std::string_view schars = text_reader::whitespace) {
		return rstrip_any(lstrip_any(s, schars), schars);
	}

	static_assert(strip_any("") == "", "strip_any");
	static_assert(strip_any(" a ") == "a", "strip_any");
	static_assert(rstrip_any(" a ") == " a", "rstrip_any");
	static_assert(lstrip_any(" a ") == "a ", "lstrip_any");

}

#endif
