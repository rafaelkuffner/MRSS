
/*
* Copyright (c) 2021 Victoria University of Wellington CMIC
* @author Ben Allen
* 
*/

#include <cstring>
#include <algorithm>

#include "iobuffer.hpp"

namespace {
	
	using namespace iob;
	using uchar = iobuffer::uchar;
	using uchar_view = iobuffer::uchar_view;

	uchar_view uchar_view_cast(std::string_view v) noexcept {
		return uchar_view{reinterpret_cast<const uchar *>(v.data()), v.size()};
	}

	std::string_view string_view_cast(uchar_view v) noexcept {
		return std::string_view{reinterpret_cast<const char *>(v.data()), v.size()};
	}

}

namespace iob {

	endian native_integer_endian() {
		constexpr auto test = uintptr_t(endian::little) | (uintptr_t(endian::big) << (CHAR_BIT * (sizeof(uintptr_t) - 1)));
		endian x{*reinterpret_cast<const uchar *>(&test)};
		assert(x == endian::little || x == endian::big);
		return x;
	}

	endian native_float_endian() {
		constexpr auto test = -1.f;
		uchar c = *reinterpret_cast<const uchar *>(&test);
		return c ? endian::big : endian::little;
	}

	std::streampos iobuffer::tell() {
		if (m_writing) {
			return tell_impl() + m_cursor - m_begin_write;
		} else {
			return tell_impl() + m_cursor - m_end_read;
		}
	}

	void iobuffer::seek(std::streamoff i, seek_origin origin) {
		// TODO optimization for seeking within buffer
		flush();
		m_begin_write = 0;
		m_cursor = 0;
		m_end_read = 0;
		m_eof = false;
		m_report_eof = false;
		// seek can be used to switch to either read or write
		// ie. user calls seek() then either get or put
		m_reading = false;
		m_writing = false;
		seek_impl(i, origin);
	}

	std::streamsize iobuffer::gavail() const noexcept {
		return m_reading ? (m_end_read - m_cursor) : 0;
	}

	std::streamsize iobuffer::pavail() const noexcept {
		return m_writing ? (capacity() - m_cursor) : 0;
	}

	uchar_view iobuffer::peek(std::streamsize n, std::streamoff base) {
		assert(base >= 0);
		if (n < 0) return {};
		read_mode();
		n += base;
		assert(m_begin_write == m_cursor);
		if (m_cursor + n <= m_end_read) {
			auto s = uchar_view{m_buf.data() + m_cursor, size_t(n)};
			s.remove_prefix(size_t(base));
			return s;
		} else if (m_eof) {
			auto s = uchar_view{m_buf.data() + m_cursor, size_t(gavail())};
			s.remove_prefix(size_t(base));
			// report eof when returning nothing with n > 0 and base offset 0
			assert(base == 0 ? (n > 0) : true);
			if (base == 0 && s.empty()) m_report_eof = true;
			return s;
		} else {
			compact();
			assert(m_cursor == 0);
			reserve(n);
			fill();
			return peek(n);
		}
	}

	uchar_view iobuffer::get(std::streamsize n) {
		auto r = peek(n);
		m_cursor += std::streamsize(r.size());
		m_begin_write = m_cursor;
		m_report_eof = std::streamsize(r.size()) < n;
		return r;
	}

	std::streamsize iobuffer::get(uchar *buf, std::streamsize n) {
		assert(buf);
		if (n < 0) return 0;
		read_mode();
		assert(m_begin_write == m_cursor);
		const auto n0 = std::min(gavail(), n);
		const auto n1 = n - n0;
		assert(n1 >= 0);
		std::memcpy(buf, m_buf.data() + m_cursor, size_t(n0));
		m_cursor += n0;
		m_begin_write = m_cursor;
		if (n1 == 0) {
			return n0;
		} else if (m_eof) {
			m_report_eof = true;
			return n0;
		} else if (n1 >= capacity()) {
			const auto n1r = get_impl(buf + n0, n1);
			if (n1r < n1) {
				m_eof = true;
				m_report_eof = true;
				m_bad = error_impl();
			}
			return n0 + n1r;
		} else {
			compact();
			assert(m_cursor == 0);
			fill();
			return n0 + get(buf + n0, n1);
		}
	}

	void iobuffer::skip_peeked(uchar_view s) {
		assert(m_reading);
		assert(m_buf.data() + m_cursor == s.data());
		assert(std::streamsize(s.size()) <= gavail());
		m_cursor += std::streamsize(s.size());
		m_begin_write = m_cursor;
	}

	void iobuffer::put(uchar_view s) {
		assert(!m_put_len && "put_begin in progress");
		const std::streamsize n(s.size());
		if (n <= 0 || m_bad) return;
		write_mode();
		assert(m_end_read == m_cursor);
		if (m_cursor + n <= capacity()) {
			std::memcpy(m_buf.data() + m_cursor, s.data(), n);
			m_cursor += n;
			m_end_read = m_cursor;
		} else {
			flush();
			compact();
			assert(m_cursor == 0 || m_bad);
			reserve(n);
			put(s);
		}
	}

	uchar * iobuffer::put_begin(std::streamsize n) {
		assert(!m_put_len && "put_begin in progress");
		if (n <= 0 || m_bad) return nullptr;
		write_mode();
		assert(m_end_read == m_cursor);
		if (pavail() >= n) {
			m_put_len = n;
			return m_buf.data() + m_cursor;
		} else {
			flush();
			compact();
			assert(m_cursor == 0 || m_bad);
			reserve(n);
			return put_begin(n);
		}
	}

	void iobuffer::put_end(std::streamsize n) {
		assert(m_put_len && "put_begin not in progress");
		assert(n >= 0 && n <= m_put_len && "overflow");
		assert(!m_bad && m_writing);
		m_cursor += n;
		m_end_read = m_cursor;
		m_put_len = 0;
	}

	void iobuffer::flush() {
		if (m_bad) return;
		if (m_begin_write < m_cursor) {
			m_begin_write += put_impl(m_buf.data() + m_begin_write, m_cursor - m_begin_write);
			if (m_begin_write < m_cursor) {
				m_bad = true;
				assert(error_impl());
			}
		}
	}

	void iobuffer::init_open() {
		m_buf.resize(default_bufsize);
	}

	void iobuffer::read_mode() {
		m_reading = true;
		if (!m_writing) return;
		// only if we were actually writing
		// flush ensures tell_impl() and m_cursor agree
		flush();
		// restart reading from cursor (because thats where tell_impl() is)
		m_end_read = m_cursor;
		m_writing = false;
	}

	void iobuffer::write_mode() {
		m_writing = true;
		if (!m_reading) return;
		// only if we were actually reading
		// seek ensures tell_impl() and m_cursor agree
		seek(tell(), seek_origin::set);
		m_reading = false;
	}

	void iobuffer::reserve(std::streamsize cap) {
		static constexpr std::streamsize blockmask = 4095;
		cap = (cap + blockmask) & blockmask;
		if (cap <= capacity()) return;
		std::vector<uchar> buf2(static_cast<size_t>(cap));
		std::memcpy(buf2.data(), m_buf.data(), size_t(m_end_read));
		m_buf = std::move(buf2);
	}

	void iobuffer::fill() {
		if (m_eof) return;
		m_end_read += get_impl(m_buf.data() + m_end_read, capacity() - m_end_read);
		if (m_end_read < capacity()) {
			m_eof = true;
			m_bad = error_impl();
		}
	}

	void iobuffer::compact() noexcept {
		std::memmove(m_buf.data(), m_buf.data() + m_begin_write, size_t(m_end_read - m_begin_write));
		m_cursor -= m_begin_write;
		m_end_read -= m_begin_write;
		m_begin_write = 0;
		// TODO shrink oversize buffer?
	}

	std::streamsize iobuffer::capacity() const noexcept {
		return std::streamsize(m_buf.size());
	}

	stream_buffer::~stream_buffer() {
		// non-owning
		m_buf = nullptr;
		m_bad = false;
	}

	stream_buffer::stream_buffer(stream_buffer &&other) noexcept :
		iobuffer(std::move(other)),
		m_buf{std::exchange(other.m_buf, nullptr)},
		m_bad{std::exchange(other.m_bad, false)}
	{}

	stream_buffer & stream_buffer::operator=(stream_buffer &&other) noexcept {
		iobuffer::operator=(std::move(other));
		m_buf = std::exchange(other.m_buf, nullptr);
		m_bad = std::exchange(other.m_bad, false);
		return *this;
	}

	stream_buffer::stream_buffer(std::streambuf *buf_) :
		m_buf{buf_}
	{
		if (m_buf) init_open();
	}

	bool stream_buffer::error_impl() {
		assert(m_buf);
		return m_bad;
	}

	void stream_buffer::seek_impl(std::streamoff i, seek_origin origin) {
		assert(m_buf);
		switch (origin) {
		case seek_origin::set:
			m_buf->pubseekpos(i);
			return;
		case seek_origin::cur:
			m_buf->pubseekoff(i, std::ios::cur);
			return;
		case seek_origin::end:
			m_buf->pubseekoff(i, std::ios::end);
			return;
		default:
			assert(false);
			return;
		}
	}

	std::streampos stream_buffer::tell_impl() {
		assert(m_buf);
		return m_buf->pubseekoff(0, std::ios::cur);
	}

	std::streamsize stream_buffer::get_impl(uchar *buf, std::streamsize n) {
		assert(m_buf);
		assert(n > 0);
		auto m = m_buf->sgetn(reinterpret_cast<char *>(buf), n);
		return m;
	}

	std::streamsize stream_buffer::put_impl(const uchar *buf, std::streamsize n) {
		assert(m_buf);
		assert(n > 0);
		auto m = m_buf->sputn(reinterpret_cast<const char *>(buf), n);
		m_buf->pubsync();
		if (m < n) m_bad = true;
		return m;
	}

	file_buffer::~file_buffer() {
		close();
	}

	file_buffer::file_buffer(file_buffer &&other) noexcept :
		iobuffer(std::move(other)),
		m_fp{std::exchange(other.m_fp, nullptr)}
	{}

	file_buffer & file_buffer::operator=(file_buffer &&other) noexcept {
		close();
		iobuffer::operator=(std::move(other));
		m_fp = std::exchange(other.m_fp, nullptr);
		return *this;
	}

	file_buffer::file_buffer(const std::filesystem::path &fpath, openmode mode) {
#ifdef _WIN32
		wchar_t zmode[8]{};
		std::wcsncpy(zmode, mode.data(), std::min<size_t>(mode.size(), 7));
		zmode[7] = L'\0';
		m_fp = _wfopen(fpath.c_str(), zmode);
#else
		char zmode[8]{};
		std::strncpy(zmode, mode.data(), std::min<size_t>(mode.size(), 7));
		zmode[7] = '\0';
		m_fp = std::fopen(fpath.c_str(), zmode);
#endif
		if (m_fp) {
			// unbuffered
			// TODO what if setvbuf fails?
			setvbuf(m_fp, nullptr, _IONBF, 0);
			init_open();
		}
	}

	void file_buffer::close() noexcept {
		if (!m_fp) return;
		try {
			flush();
		} catch (...) {}
		std::fclose(m_fp);
		m_fp = nullptr;
	}

	bool file_buffer::error_impl() {
		assert(m_fp);
		return bool(std::ferror(m_fp));
	}

	void file_buffer::seek_impl(std::streamoff i, seek_origin origin) {
		assert(m_fp);
		// TODO what if fseek fails?
#ifdef _WIN32
		_fseeki64(m_fp, i, int(origin));
#else
		static_assert(sizeof(std::streamoff) == sizeof(long), "fseek: sizeof(offset)");
		std::fseek(m_fp, i, int(origin));
#endif
	}

	std::streampos file_buffer::tell_impl() {
		assert(m_fp);
#ifdef _WIN32
		return _ftelli64(m_fp);
#else
		return std::ftell(m_fp);
#endif
	}

	std::streamsize file_buffer::get_impl(uchar *buf, std::streamsize n) {
		assert(m_fp);
		assert(n > 0);
		return std::streamoff(std::fread(buf, 1, size_t(n), m_fp));
	}

	std::streamsize file_buffer::put_impl(const uchar *buf, std::streamsize n) {
		assert(m_fp);
		assert(n > 0);
		auto r = std::streamsize(std::fwrite(buf, 1, size_t(n), m_fp));
		std::fflush(m_fp);
		return r;
	}

	void text_reader::seek(std::streamoff i, iobuffer::seek_origin origin) {
		assert(m_iobuf);
		m_iobuf->seek(i, origin);
	}

	std::streampos text_reader::tell() const {
		assert(m_iobuf);
		return m_iobuf->tell();
	}

	std::string_view text_reader::peek(std::streamsize n, std::streamoff base) {
		assert(m_iobuf);
		auto v = m_iobuf->peek(n, base);
		return string_view_cast(v);
	}

	std::string_view text_reader::get(std::streamsize n) {
		assert(m_iobuf);
		auto v = m_iobuf->get(n);
		return string_view_cast(v);
	}

	void text_reader::skip_peeked(std::string_view s) {
		assert(m_iobuf);
		m_iobuf->skip_peeked(uchar_view_cast(s));
	}

	std::string_view text_reader::peek_until_any(std::string_view delimchars, std::streamoff base) {
		// TODO allow peek whole input?
		assert(delimchars.size());
		return peek_until([=](std::string_view v, std::string::size_type i) {
			return v.find_first_of(delimchars, i);
		}, base);
	}

	std::string_view text_reader::peek_until_any(char delim, std::streamoff base) {
		return peek_until([=](std::string_view v, std::string::size_type i) {
			return v.find(delim, i);
		}, base);
	}

	std::string_view text_reader::peek_until_seq(std::string_view delim, std::streamoff base) {
		return peek_until([=](std::string_view v, std::string::size_type i) {
			// note: need to check overlap of delim with old/new chars
			// TODO test this
			return v.find(delim, std::max(i + 1, delim.size()) - delim.size());
		}, base);
	}

	std::string_view text_reader::peek_while_any(std::string_view tokchars, std::streamoff base) {
		if (!tokchars.size()) return {};
		return peek_until([=](std::string_view v, std::string::size_type i) {
			return v.find_first_not_of(tokchars, i);
		}, base);
	}

	std::string_view text_reader::peek_while_any(char tokchar, std::streamoff base) {
		return peek_until([=](std::string_view v, std::string::size_type i) {
			return v.find_first_not_of(tokchar, i);
		}, base);
	}

	std::string_view text_reader::get_until_any(std::string_view delimchars) {
		auto s = peek_until_any(delimchars);
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::get_until_any(char delim) {
		auto s = peek_until_any(delim);
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::get_until_seq(std::string_view delim) {
		auto s = peek_until_seq(delim);
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::get_while_any(std::string_view tokchars) {
		auto s = peek_while_any(tokchars);
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::get_while_any(char tokchar) {
		auto s = peek_while_any(tokchar);
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::get_line(char delim, char strip) {
		// peek with 1 extra char to include the delim
		auto s = peek_until([=](std::string_view v, std::string::size_type i) {
			return v.find(delim, i);
		}, 0, 1);
		skip_peeked(s);
		if (s.size() && s.back() == delim) s.remove_suffix(1);
		if (s.size() && s.back() == strip) s.remove_suffix(1);
		return s;
	}

	void text_reader::skip_line(char delim) {
		// peek with 1 extra char to include the delim
		auto s = peek_until([=](std::string_view v, std::string::size_type i) {
			return v.find(delim, i);
		}, 0, 1);
		skip_peeked(s);
	}

	bool text_reader::has_any(std::string_view tokchars) {
		auto s = peek(1);
		if (s.empty()) return false;
		return tokchars.find(s[0]) != std::string::npos;
	}

	bool text_reader::has_any(char c) {
		auto s = peek(1);
		if (s.empty()) return false;
		return s[0] == c;
	}

	bool text_reader::has_seq(std::string_view tok) {
		auto s = peek(tok.size());
		return s == tok;
	}

	std::string_view text_reader::get_any(std::string_view tokchars) {
		auto s = peek(1);
		if (s.empty()) return {};
		if (tokchars.find(s[0]) == std::string::npos) return {};
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::get_any(char c) {
		auto s = peek(1);
		if (s.empty()) return {};
		if (s[0] != c) return {};
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::get_seq(std::string_view tok) {
		auto s = peek(tok.size());
		if (s != tok) return {};
		skip_peeked(s);
		return s;
	}

	std::string_view text_reader::peek_until_ws() {
		return peek_until_any(whitespace);
	}

	std::string_view text_reader::peek_while_ws() {
		return peek_while_any(whitespace);
	}

	std::string_view text_reader::get_until_ws() {
		return get_until_any(whitespace);
	}

	std::string_view text_reader::get_while_ws() {
		return get_while_any(whitespace);
	}

	void text_reader::skip_ws() {
		get_while_ws();
	}

	void binary_reader::seek(std::streamoff i, iobuffer::seek_origin origin) {
		assert(m_iobuf);
		m_iobuf->seek(i, origin);
	}

	std::streampos binary_reader::tell() const {
		assert(m_iobuf);
		return m_iobuf->tell();
	}

	uchar_view binary_reader::peek(std::streamsize n, std::streamoff base) {
		assert(m_iobuf);
		return m_iobuf->peek(n, base);
	}

	uchar_view binary_reader::get(std::streamsize n) {
		assert(m_iobuf);
		return m_iobuf->get(n);
	}

	void binary_reader::skip_peeked(uchar_view s) {
		assert(m_iobuf);
		m_iobuf->skip_peeked(s);
	}

	std::string_view binary_reader::get_str(std::streamsize n) {
		return string_view_cast(get(n));
	}

	void text_writer::seek(std::streamoff i, iobuffer::seek_origin origin) {
		assert(m_iobuf);
		m_iobuf->seek(i, origin);
	}

	std::streampos text_writer::tell() const {
		assert(m_iobuf);
		return m_iobuf->tell();
	}

	void text_writer::flush() {
		assert(m_iobuf);
		m_iobuf->flush();
	}

	void text_writer::put(char c) {
		char *p = put_begin(1);
		*p = c;
		put_end(1);
	}

	void text_writer::put(std::string_view s) {
		assert(m_iobuf);
		m_iobuf->put(uchar_view_cast(s));
	}

	char * text_writer::put_begin(std::streamsize n) {
		assert(m_iobuf);
		return reinterpret_cast<char *>(m_iobuf->put_begin(n));
	}

	void text_writer::put_end(std::streamsize n) {
		assert(m_iobuf);
		m_iobuf->put_end(n);
	}

	void binary_writer::seek(std::streamoff i, iobuffer::seek_origin origin) {
		assert(m_iobuf);
		m_iobuf->seek(i, origin);
	}

	std::streampos binary_writer::tell() const {
		assert(m_iobuf);
		return m_iobuf->tell();
	}

	void binary_writer::flush() {
		assert(m_iobuf);
		m_iobuf->flush();
	}

	void binary_writer::put(uchar c) {
		uchar *p = put_begin(1);
		*p = c;
		put_end(1);
	}

	void binary_writer::put(uchar_view s) {
		assert(m_iobuf);
		m_iobuf->put(s);
	}

	uchar * binary_writer::put_begin(std::streamsize n) {
		assert(m_iobuf);
		return m_iobuf->put_begin(n);
	}

	void binary_writer::put_end(std::streamsize n) {
		assert(m_iobuf);
		m_iobuf->put_end(n);
	}

	void binary_writer::put_str(std::string_view s) {
		assert(m_iobuf);
		m_iobuf->put(uchar_view_cast(s));
	}

}
