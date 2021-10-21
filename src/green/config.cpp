
#include "config.hpp"

#include "iobuffer/iobuffer.hpp"

#include <iostream>

namespace green {

	bool load_config(const std::filesystem::path &fpath) {
		auto &presets = saliency_presets();
		iob::file_buffer fbuf(fpath, iob::file_buffer::read);
		if (!fbuf.is_open()) {
			std::cout << "failed to open " << fpath.u8string() << std::endl;
			return false;
		}
		std::cout << "loading saliency presets from " << fpath.u8string() << std::endl;
		iob::text_reader in(&fbuf);
		while (in.good()) {
			in.skip_ws();
			auto funcname = in.get_until_ws();
			in.skip_ws();
			if (funcname.empty() || funcname[0] == '#') {
				in.skip_line();
			} else if (funcname == "sal_preset") {
				saliency_preset p;
				p.name = in.get_until_ws();
				in.skip_ws();
				auto paramstr = in.get_until_ws();
				if (p.uparams.parse(paramstr)) {
					std::cout << "parsed saliency preset " << p.name << " : " << p.uparams.str() << std::endl;
					presets.push_back(std::move(p));
				} else {
					std::cout << "failed to parse saliency preset " << p.name << std::endl;
				}
				in.skip_line();
			} else {
				std::cout << "unknown config function " << funcname << std::endl;
				in.skip_line();
			}
		}
		return true;
	}

	bool save_saliency_presets(iob::iobuffer &buf, const std::vector<saliency_preset> &presets) {
		iob::text_writer out(&buf);
		out.put("# MRSS saliency presets\n");
		for (auto &p : presets) {
			if (p.builtin) continue;
			out.put("sal_preset ");
			out.put(p.name);
			out.put(' ');
			out.put(p.uparams.str());
			out.put('\n');
		}
		return out.good();
	}

	bool save_presets(const std::filesystem::path &fpath) {
		iob::file_buffer fbuf(fpath, iob::file_buffer::write);
		if (!fbuf.is_open()) {
			std::cout << "failed to open " << fpath.u8string() << std::endl;
			return false;
		}
		std::cout << "saving presets to " << fpath.u8string() << std::endl;
		bool r = true;
		r |= save_saliency_presets(fbuf, saliency_presets());
		return r;
	}

	bool save_saliency_presets(const std::filesystem::path &fpath, const std::vector<saliency_preset> &presets) {
		iob::file_buffer fbuf(fpath, iob::file_buffer::write);
		if (!fbuf.is_open()) {
			std::cout << "failed to open " << fpath.u8string() << std::endl;
			return false;
		}
		std::cout << "saving saliency presets to " << fpath.u8string() << std::endl;
		return save_saliency_presets(fbuf, presets);
	}

}
