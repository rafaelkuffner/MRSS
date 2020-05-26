
#ifdef GREEN_DIALOG_WIN32

#define NOMINMAX
#include <Windows.h>

#include <thread>
#include <algorithm>
#include <iostream>

#include "dialog.hpp"

#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>

namespace stdfs = std::filesystem;

namespace green {

	std::future<std::vector<std::filesystem::path>> open_file_dialog(GLFWwindow *win, const std::filesystem::path &hint, bool multi) {
		const HWND hwnd = glfwGetWin32Window(win);
		return std::async([=]() {
			constexpr int nbuf = 1024;
			wchar_t dirbuf[nbuf]{};
			wchar_t buf[nbuf]{};
			if (stdfs::is_directory(hint)) {
				lstrcpynW(dirbuf, hint.c_str(), nbuf);
			} else {
				lstrcpynW(buf, hint.c_str(), nbuf);
			}
			OPENFILENAMEW ofn{};
			ofn.lStructSize = sizeof(OPENFILENAMEW);
			ofn.hwndOwner = hwnd;
			ofn.lpstrInitialDir = dirbuf;
			ofn.lpstrFile = buf;
			ofn.nMaxFile = nbuf;
			ofn.Flags |= OFN_EXPLORER | OFN_FILEMUSTEXIST;
			if (multi) ofn.Flags |= OFN_ALLOWMULTISELECT;
			if (!GetOpenFileNameW(&ofn)) {
				const DWORD err = CommDlgExtendedError();
				if (err != 0) {
					std::cerr << "win32 dialog error " << err << std::endl;
				}
				return std::vector<stdfs::path>{};
			}
			std::vector<stdfs::path> v;
			stdfs::path dir;
			// parse file paths - null separated, first is dirpath, then filenames
			// if only one filename then only one element, no null separation
			for (wchar_t *p0 = buf, *p1 = buf; *p0 && (p1 = std::find(p0, buf + nbuf, L'\0')) < buf + nbuf - 1; p0 = p1 + 1) {
				if (p0 == buf && *(p1 + 1)) {
					// first elem - dir
					dir.assign(p0, p1);
				} else if (p0 == buf) {
					// first elem, no separator - single file
					v.emplace_back(p0, p1);
				} else {
					// filename only
					v.push_back(dir / stdfs::path(p0, p1));
				}
			}
			std::cerr << "open files:" << std::endl;
			for (const auto &p : v) {
				std::cerr << p.u8string() << std::endl;
			}
			return v;
		});
	}

	std::future<std::filesystem::path> save_file_dialog(GLFWwindow *win, const std::filesystem::path &hint) {
		const HWND hwnd = glfwGetWin32Window(win);
		return std::async([=]() {
			constexpr int nbuf = 1024;
			wchar_t dirbuf[nbuf]{};
			wchar_t buf[nbuf]{};
			if (stdfs::is_directory(hint)) {
				lstrcpynW(dirbuf, hint.c_str(), nbuf);
			} else {
				lstrcpynW(buf, hint.c_str(), nbuf);
			}
			OPENFILENAMEW ofn{};
			ofn.lStructSize = sizeof(OPENFILENAMEW);
			ofn.hwndOwner = hwnd;
			ofn.lpstrInitialDir = dirbuf;
			ofn.lpstrFile = buf;
			ofn.nMaxFile = nbuf;
			ofn.Flags |= OFN_EXPLORER | OFN_OVERWRITEPROMPT;
			if (!GetSaveFileNameW(&ofn)) {
				const DWORD err = CommDlgExtendedError();
				if (err != 0) {
					std::cerr << "win32 dialog error " << err << std::endl;
				}
				return stdfs::path();
			}
			stdfs::path p{buf};
			std::cerr << "save file:" << p.u8string() << std::endl;
			return p;
		});
	}

}

#endif
