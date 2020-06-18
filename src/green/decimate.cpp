
/*
* Copyright 2020
* Computational Media Innovation Centre
* Victoria University of Wellington
*
*/

#include <iostream>

#include "OpenMesh/Tools/Decimater/DecimaterT.hh"
#include "OpenMesh/Tools/Decimater/ModQuadricT.hh"
#include "ModBoundedSaliencyT.hh"
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

	

		lce::ModBoundedSaliencyT<TriMesh>::Handle hModBounded;
		decimater.add(hModBounded);
		decimater.module(hModBounded).setSaliencyProperty(mparams.prop_saliency);
		
		int step1RemovedVerts = (mparams.mesh->n_vertices() - uparams.targetverts) * 0.6;
		int step2RemovedVerts = (mparams.mesh->n_vertices() - uparams.targetverts) * 0.3;
		int step3RemovedVerts = (mparams.mesh->n_vertices() - uparams.targetverts) * 0.1;



		// TODO how do i get some sort of progress indicator?
		decimater.initialize();

		decimater.module(hModBounded).setMinCollapsePercentile(0);
		decimater.module(hModBounded).setMaxCollapsePercentile(0.33);
		std::cout << "Decimating 1st step to " << mparams.mesh->n_vertices() - step1RemovedVerts << " vertices" << std::endl;
		decimater.decimate_to(mparams.mesh->n_vertices() - step1RemovedVerts);
		mparams.mesh->garbage_collection();

		decimater.module(hModBounded).setMinCollapsePercentile(0.34);
		decimater.module(hModBounded).setMaxCollapsePercentile(0.66); 
		std::cout << "Decimating 2nd step to " << mparams.mesh->n_vertices() - step2RemovedVerts << " vertices" << std::endl;
		decimater.initialize();
		decimater.decimate_to(mparams.mesh->n_vertices() - step2RemovedVerts);
		mparams.mesh->garbage_collection();

		decimater.module(hModBounded).setMinCollapsePercentile(0.67);
		decimater.module(hModBounded).setMaxCollapsePercentile(1.00);
		std::cout << "Decimating 3rd step to " << mparams.mesh->n_vertices() - step3RemovedVerts << " vertices" << std::endl;
		decimater.initialize();
		decimater.decimate_to(mparams.mesh->n_vertices() - step3RemovedVerts);
		mparams.mesh->garbage_collection();


		progress.state = decimation_state::done;
		progress.elapsed_time = std::chrono::duration_cast<decltype(decimate_progress::elapsed_time)>(std::chrono::steady_clock::now() - time_start);
		return true;
	}

}
