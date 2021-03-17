
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
#include <vector>
#include <tuple>
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
			return m_cursor == m_end_read && m_eof;
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
		std::string_view peek_until(FindDelim &&find, std::streamoff base = 0) {
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
						return peek(i, base);
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

		std::string_view get_until_ws();

		std::string_view get_while_ws();

		template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
		IOB_NOINLINE
		std::errc get_int(T &val, int base = 10) {
			auto s = peek(32);
			auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), val, base);
			get(end - s.data());
			return ec;
		}

		template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
		IOB_NOINLINE
		std::errc get_float(T &val, std::chars_format fmt = std::chars_format::general) {
			auto s = peek(32);
			auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), val, fmt);
			get(end - s.data());
			return ec;
		}

		template <typename ...Ts, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_ints(std::tuple<Ts &...> vals, Delim delimchars = whitespace, int base = 10) {
			return get_ints_impl(vals, std::index_sequence_for<Ts...>(), delimchars, base);
		}

		template <typename ...Ts, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_floats(std::tuple<Ts &...> vals, Delim delimchars = whitespace, std::chars_format fmt = std::chars_format::general) {
			return get_floats_impl(vals, std::index_sequence_for<Ts...>(), delimchars, fmt);
		}

	private:
		template <typename ...Ts, size_t ...Is, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_ints_impl(std::tuple<Ts &...> vals, std::index_sequence<Is...>, Delim delimchars, int base = 10) {
			std::errc ec{};
			int r = 0;
			(... && (get_while_any(delimchars), (ec = get_int(std::get<Is>(vals), base)), r += int(ec == std::errc{}), ec == std::errc{}));
			return {ec, r};
		}

		template <typename ...Ts, size_t ...Is, typename Delim>
		IOB_FORCEINLINE
		std::pair<std::errc, int> get_floats_impl(std::tuple<Ts &...> vals, std::index_sequence<Is...>, Delim delimchars, std::chars_format fmt) {
			std::errc ec{};
			int r = 0;
			(... && (get_while_any(delimchars), (ec = get_float(std::get<Is>(vals), fmt)), r += int(ec == std::errc{}), ec == std::errc{}));
			return {ec, r};
		}

	};

}

#endif
