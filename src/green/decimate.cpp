
/*
* Copyright 2020
* Computational Media Innovation Centre
* Victoria University of Wellington
*
*/

#include <iostream>

#include "OpenMesh/Tools/Decimater/DecimaterT.hh"
#include "OpenMesh/Tools/Decimater/ModQuadricT.hh"

#include "decimate.hpp"

namespace green {
	
	void decimate_user_params::sanitize() {
		targetverts = std::max(targetverts, 0);
	}

	bool decimate(const decimate_mesh_params &mparams, const decimate_user_params &uparams, decimate_progress &progress) {
		progress.state = decimation_state::run;
		const auto time_start = std::chrono::steady_clock::now();
		
		OpenMesh::Decimater::DecimaterT<TriMesh> decimater(*mparams.mesh);
		OpenMesh::Decimater::ModQuadricT<TriMesh>::Handle hModQuadrics;
		decimater.add(hModQuadrics);
		decimater.module(hModQuadrics).unset_max_err();

		// TODO how do i get some sort of progress indicator?
		std::cout << "Decimating to " << uparams.targetverts << " vertices" << std::endl;
		decimater.initialize();
		decimater.decimate_to(uparams.targetverts);
		mparams.mesh->garbage_collection();

		progress.state = decimation_state::done;
		progress.elapsed_time = std::chrono::duration_cast<decltype(decimate_progress::elapsed_time)>(std::chrono::steady_clock::now() - time_start);
		return true;
	}

}
