
#pragma once

#include "uilocale.hpp"

namespace green {

	constexpr uilocale uilocale_init_en() {
		using namespace uistrings;
		uilocale loc;
		// common
		loc[help_input] = "Input file";
		loc[help_output] = "Output file";
		loc[help_threads] = "Number of threads to use for parallel sections.\nDefault is all available processors.";
		loc[help_cli_opts_common] = "Common options:";
		loc[help_cli_version] = "Show version";
		loc[help_cli_help] = "Show help";
		loc[help_cli_ascii] = "Save in plaintext instead of binary if possible";
		loc[help_cli_gui] = "Show result in gui when done";
		loc[help_cli_noprogress] = "Suppress progress output (useful for automation)";
		loc[param_threads] = "Threads";
		// saliency
		loc[help_cli_opts_saliency] = "Saliency options:";
		loc[help_sal_go] = "Compute saliency";
		loc[help_sal_area] = "Size of the largest salient features.\nSpecified as a fraction of the surface area.";
		loc[help_sal_levels] = "Using more levels allows increasingly local features to become visible.";
		loc[help_sal_normpower] = "Controls how quickly saliency tends towards extreme values.";
		loc[help_sal_curvweight] = "Additional visibility for immediate local features of high curvature.";
		loc[help_sal_noisefilter] = "Filter out noise from otherwise smooth surfaces.\nEnable to access noise parameters.";
		loc[help_sal_noiseheight] = "Maximum magnitude of noise to filter.\nSpecified as a fraction of the square root of the surface area.";
		loc[help_sal_samplespern] = "Higher: more samples, more accurate results.\nDoes not usually need tuning per model.";
		loc[help_sal_subsample] = "Apply automatic subsampling (fast, recommended).\nEnable to access sampling parameters.";
		loc[help_sal_full] = "Disable subsampling (slow, not recommended)";
		loc[param_sal_area] = "Area";
		loc[param_sal_levels] = "Levels";
		loc[param_sal_normpower] = "Contrast";
		loc[param_sal_curvweight] = "Contour";
		loc[param_sal_noisefilter] = "Noise Filter";
		loc[param_sal_noiseheight] = "Noise Height";
		loc[param_sal_samplespern] = "S/N";
		loc[param_sal_subsample] = "Subsample";
		// decimation
		loc[help_cli_opts_decimate] = "Decimation options (Work in Progress):";
		loc[help_dec_go] = "Perform decimation";
		loc[help_dec_targetverts] = "Target number of vertices";
		loc[help_dec_weight] = "Saliency weighting; 0 (even) .. 1 (most weight to highest saliency)";
		loc[help_dec_power] = "Non-linearity of saliency weighting; 1 (linear) .. 2 (more extreme)";
		loc[help_dec_bins] = "Number of saliency bins";
		loc[help_dec_propname] = "Mesh vertex property to use for decimation if not calculating saliency.\nDefault is 'quality'.";
		loc[help_dec_usesaliency] = "Use saliency for decimation";
		loc[help_dec_nosaliency] = "Don't use saliency for decimation";
		loc[param_dec_usesaliency] = "Use Saliency";
		loc[param_dec_targetverts] = "Vertices";
		loc[param_dec_weight] = "Weight";
		loc[param_dec_power] = "Power";
		loc[param_dec_bins] = "Bins";
		return loc;
	}

	const uilocale & uilocale_en() {
		static constexpr uilocale loc = uilocale_init_en();
		static_assert(!loc.any_empty(), "empty string in localization");
		return loc;
	}

}
