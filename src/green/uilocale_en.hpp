
#pragma once

#include "uilocale.hpp"

namespace green {

	constexpr uilocale uilocale_init_en() {
		using namespace uistrings;
		uilocale loc;
		// common
		loc[help_input] = "Input file";
		loc[help_output] = "Output file";
		loc[help_threads] = "Number of threads to use for parallel sections.\nDefault (0) is all available processors.";
		loc[help_cli_opts_common] = "Common options:";
		loc[help_cli_version] = "Show version";
		loc[help_cli_help] = "Show help";
		loc[help_cli_ascii] = "Save in plaintext instead of binary if possible";
		loc[help_cli_original_vids] = "Save extra vertex property with original pre-decimation vertex ids";
		loc[help_cli_gui] = "Show result in gui when done";
		loc[help_cli_noprogress] = "Suppress progress output (useful for automation)";
		loc[help_cli_color] = "How to colorize the model (vertex colors) when saving. Options are:\n"
			"- none (don't colorize)\n- vcolor (preserve original vertex colors)\n- saliency (colorize a saliency property)\n";
			//"- saliency_comparison (colorize difference of saliency properties)\n";
		loc[help_cli_propnames] = "    Properties can be specified by name (as seen in the UI; no 'quality' prefix or metadata)\n"
			"    or by index with @index (starting at @0). Negative indices can be used for reverse indexing from\n"
			"    the last property (which is @-1). @index notation cannot be used when assigning a property name\n"
			"    (as with --salprop); property indices are automatically assigned when saving the model.";
		loc[help_cli_config] = "Config file to load for e.g. saliency presets";
		loc[param_threads] = "Threads";
		// saliency
		loc[help_cli_opts_saliency] = "Saliency options:";
		loc[help_sal_go] = "Compute saliency";
		loc[help_sal_curv] = "Curvature measure to use. Options are:\n"
			"- don (default)\n- mean (for comparison with previous work, disabled by default)";
		loc[help_sal_area] = "Size of the largest salient features.\nSpecified as a fraction of the surface area.";
		loc[help_sal_levels] = "Using more levels allows increasingly local features to become visible.";
		loc[help_sal_normpower] = "Controls how quickly saliency tends towards extreme values.\n"
			"Low/high values cause changes around low/high curvatures to be more salient.";
		loc[help_sal_autocontrast] = "Automatically set contrast to make mean saliency approximately midrange (white).";
		loc[help_sal_curvweight] = "Additional visibility for immediate local features of high curvature.";
		loc[help_sal_noisefilter] = "Filter out noise from otherwise smooth surfaces.\nEnable to access noise parameters.";
		loc[help_sal_noiseheight] = "Maximum magnitude of noise to filter.\nSpecified as a fraction of the square root of the surface area.";
		loc[help_sal_samplespern] = "Higher: more samples, more accurate results.\nDoes not usually need tuning per model.";
		loc[help_sal_subsample] = "Apply automatic subsampling (fast, recommended).\nEnable to access sampling parameters.";
		loc[help_sal_preset] = "Saliency preset to use";
		loc[help_sal_full] = "Disable subsampling (slow, not recommended)";
		loc[help_sal_salprop] = "Name for computed saliency property.\n"
			"Will always be placed at the last index (@-1) for subsequent decimation, colorization and saving.";
		loc[help_sal_colorprop] = "Saliency property to colorize when saving.\n"
			"Default is the saliency property used for decimation (--decprop), or the computed saliency (@-1), or @0.";
		loc[help_sal_globality] = "Single meta parameter to adjust saliency between local and global.";
		loc[param_sal_curv] = "Curvature";
		loc[param_sal_area] = "Area";
		loc[param_sal_levels] = "Levels";
		loc[param_sal_normpower] = "Contrast";
		loc[param_sal_curvweight] = "Contour";
		loc[param_sal_noisefilter] = "Noise Filter";
		loc[param_sal_noiseheight] = "Noise Height";
		loc[param_sal_samplespern] = "S/N";
		loc[param_sal_subsample] = "Subsample";
		loc[param_sal_autocontrast] = "Automatic contrast";
		loc[param_sal_globality] = "Globality";
		// decimation
		loc[help_cli_opts_decimate] = "Decimation options (Work in Progress):";
		loc[help_dec_go] = "Perform decimation";
		loc[help_dec_usebins] = "Use saliency bins instead of weighting";
		loc[help_dec_usetris] = "Decimate by triangle count (not exact) instead of vertex count (exact)";
		loc[help_dec_targetverts] = "Target number of vertices (exact)";
		loc[help_dec_targettris] = "Target number of triangles (not exact)";
		loc[help_dec_weight] = "Saliency weighting (error multiplier):\n"
			"-  k : vertices with sal=1 will have 1/k error compared to those with sal=0\n"
			"-  1 : even weighting, saliency has no impact\n"
			"- <1 : more salient => more error\n"
			"- >1 : more salient => less error\n"
			"- 20 : default"; // TODO keep default value up-to-date
		loc[help_dec_binweight] = "Saliency bin weighting; 0 (even) .. 1 (most weight to highest saliency)";
		loc[help_dec_power] = "Non-linearity of saliency weighting:\n"
			"- <1 : mid-range saliency preserved more\n"
			"-  1 : linear\n"
			"- >1 : mid-range saliency preserved less\n"
			"-  2 : default"; // TODO keep default value up-to-date
		loc[help_dec_bins] = "Number of saliency bins";
		loc[help_dec_seams] = "Preserve seams in UV coords and normals";
		loc[help_dec_folds] = "Prevent mesh surface folding over itself (flipped triangles)";
		loc[help_dec_limit_aspect] = "Limit triangle aspect ratio";
		loc[help_dec_max_aspect] = "Maximum desired triangle aspect ratio";
		loc[help_dec_decprop] = "Saliency property to use for decimation if not calculating saliency.\n"
			"Default is @0, the first property.";
		loc[help_dec_usesaliency] = "Use saliency for decimation";
		loc[help_dec_nosaliency] = "Don't use saliency for decimation";
		loc[param_dec_usesaliency] = "Use Saliency";
		loc[param_dec_usebins] = "Use Bins";
		loc[param_dec_usetris] = "Use Triangles";
		loc[param_dec_targetverts] = "Vertices";
		loc[param_dec_targettris] = "Triangles";
		loc[param_dec_weight] = "Sal Weight";
		loc[param_dec_binweight] = "Bin Weight";
		loc[param_dec_power] = "Sal Power";
		loc[param_dec_binpower] = "Bin Power";
		loc[param_dec_bins] = "Bins";
		loc[param_dec_seams] = "Preserve Seams";
		loc[param_dec_folds] = "Prevent Folds";
		loc[param_dec_limit_aspect] = "Limit Aspect Ratio";
		loc[param_dec_max_aspect] = "Max Aspect";
		// curvature
		loc[curv_don] = "DoN curvature";
		loc[curv_mean] = "Mean curvature";
		return loc;
	}

	const uilocale & uilocale_en() {
		static constexpr uilocale loc = uilocale_init_en();
		static_assert(!loc.any_empty(), "empty string in localization");
		return loc;
	}

}
