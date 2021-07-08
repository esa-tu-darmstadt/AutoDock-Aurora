#include "auxiliary.h"

void energy_and_gradient (
	const	uint 				DockConst_pop_size,
			float*				final_interE,
			float*				final_intraE,

			float*				local_coords_x,
			float*				local_coords_y,
			float* 				local_coords_z,

/*
			float				local_coords_x[][MAX_POPSIZE],
			float				local_coords_y[][MAX_POPSIZE],
			float 				local_coords_z[][MAX_POPSIZE],
*/
			float*				gradient_inter_x,
			float*				gradient_inter_y,
			float*				gradient_inter_z,
			float*				gradient_intra_x,
			float*				gradient_intra_y,
			float*				gradient_intra_z,
	// ie
	const 	float* 	restrict	IE_Fgrids,
	const 	float*	restrict	IA_IE_atom_charges,
	const	int*	restrict	IA_IE_atom_types,
			uchar 				DockConst_xsz,
			uchar 				DockConst_ysz,
			uchar 				DockConst_zsz,
			uchar 				DockConst_num_of_atoms,	// TODO: FIXME
			float 				DockConst_gridsize_x_minus1,
			float 				DockConst_gridsize_y_minus1,
			float 				DockConst_gridsize_z_minus1,
			uint 				Host_mul_tmp2,
			uint 				Host_mul_tmp3,
	// ia
	const 	int*	restrict	IA_intraE_contributors,
	const	float*	restrict	IA_reqm,
	const 	float*	restrict	IA_reqm_hbond,
	const 	uint*	restrict	IA_atom1_types_reqm,
	const	uint*	restrict	IA_atom2_types_reqm,
	const	float*	restrict	IA_VWpars_AC,
	const	float*	restrict	IA_VWpars_BD,
	const	float*	restrict	IA_dspars_S,
	const 	float*	restrict	IA_dspars_V,
			float				DockConst_smooth,
			uint 				DockConst_num_of_intraE_contributors,
			float 				DockConst_grid_spacing,
			uint 				DockConst_num_of_atypes,
			float				DockConst_coeff_elec,
			float				DockConst_qasp,
			float				DockConst_coeff_desolv,
    // gradients
    const   int*    restrict    GRAD_rotbonds,
    const   int*    restrict    GRAD_rotbonds_atoms,
    const   int*    restrict    GRAD_num_rotating_atoms_per_rotbond,

    const   float*  restrict    GRAD_angle,
    const   float*  restrict    GRAD_dependence_on_theta,
    const   float*  restrict    GRAD_dependence_on_rotangle
) {
	// Interpreting IE_Fgrids as multidimensional static array
	const float (*IE_Fg)[MAX_NUM_OF_ATYPES+2][DockConst_zsz][DockConst_ysz][DockConst_xsz] =
		(void*)(IE_Fgrids);
	const float (*IE_Fg_2)[DockConst_zsz][DockConst_ysz][DockConst_xsz] =
		(void*)(&IE_Fgrids[Host_mul_tmp2]);
	const float (*IE_Fg_3)[DockConst_zsz][DockConst_ysz][DockConst_xsz] =
		(void*)(&IE_Fgrids[Host_mul_tmp3]);


	for (uint atom_id = 0; atom_id < DockConst_num_of_atoms; atom_id++) {

		int atom_typeid = IA_IE_atom_types[atom_id];
		float q = IA_IE_atom_charges[atom_id];

		float x = local_coords_x[atom_id];
		float y = local_coords_y[atom_id];
		float z = local_coords_z[atom_id];

		float partialE1;
		float partialE2;
		float partialE3;

		// If atom is outside of the grid
		if ((x < 0.0f) || (x >= DockConst_gridsize_x_minus1) ||
			(y < 0.0f) || (y >= DockConst_gridsize_y_minus1) ||
			(z < 0.0f) || (z >= DockConst_gridsize_z_minus1)) {

			x = x - 0.5f * DockConst_gridsize_x_minus1; // TODO: in ADGPU is just DockConst_gridsize_x
			y = y - 0.5f * DockConst_gridsize_y_minus1;
			z = z - 0.5f * DockConst_gridsize_z_minus1;

			partialE1 = 16777216.0f;
			partialE2 = 0.0f;
			partialE3 + 0.0f;

			gradient_inter_x[atom_id] = 16777216.0f;
			gradient_inter_y[atom_id] = 16777216.0f;
			gradient_inter_z[atom_id] = 16777216.0f;
		}
		else {
			float x_low = floorf(x);
			float y_low = floorf(y);
			float z_low = floorf(z);

			int ix = (int)x_low;
			int iy = (int)y_low;
			int iz = (int)z_low;

			float dx = x - x_low;
			float dy = y - y_low;
			float dz = z - z_low;

			float omdx = 1.0f - dx;
			float omdy = 1.0f - dy;
			float omdz = 1.0f - dz;

			// Calculating interpolation weights
			float weight000, weight001, weight010;
			float weight011, weight100, weight101;
			float weight110, weight111;
			weight000 = omdx * omdy * omdz;
			weight100 = dx * omdy * omdz;
			weight010 = omdx * dy * omdz;
			weight110 = dx * dy * omdz;
			weight001 = omdx * omdy * dz;
			weight101 = dx * omdy * dz;
			weight011 = omdx * dy * dz;
			weight111 = dx * dy * dz;

#ifdef (PRINT_ALL)
			printf("\n\nPartial results for atom with id %i:\n", atom_id);
			printf("x_low = %f, x_frac = %f\n", x_low, dx);
			printf("y_low = %f, y_frac = %f\n", y_low, dy);
			printf("z_low = %f, z_frac = %f\n", z_low, dz);
			printf("weight(0,0,0) = %f\n", weight000);
			printf("weight(1,0,0) = %f\n", weight100);
			printf("weight(0,1,0) = %f\n", weight010);
			printf("weight(1,1,0) = %f\n", weight110);
			printf("weight(0,0,1) = %f\n", weight001);
			printf("weight(1,0,1) = %f\n", weight101);
			printf("weight(0,1,1) = %f\n", weight011);
			printf("weight(1,1,1) = %f\n", weight111);
#endif

			// Energy contribution of the current grid type
			float cub000, cub001, cub010;
			float cub011, cub100, cub101;
			float cub110, cub111;

			cub000 = (*IE_Fg)[atom_typeid][iz  ][iy  ][ix  ];
			cub100 = (*IE_Fg)[atom_typeid][iz  ][iy  ][ix+1];
			cub010 = (*IE_Fg)[atom_typeid][iz  ][iy+1][ix  ];
			cub110 = (*IE_Fg)[atom_typeid][iz  ][iy+1][ix+1];
			cub001 = (*IE_Fg)[atom_typeid][iz+1][iy  ][ix  ];
			cub100 = (*IE_Fg)[atom_typeid][iz+1][iy  ][ix+1];
			cub011 = (*IE_Fg)[atom_typeid][iz+1][iy+1][ix  ];
			cub111 = (*IE_Fg)[atom_typeid][iz+1][iy+1][ix+1];

#ifdef (PRINT_ALL)
			printf("Interpolation of van der Waals map:\n");
			printf("cube(0,0,0) = %f\n". cub000);
			printf("cube(1,0,0) = %f\n". cub100);
			printf("cube(0,1,0) = %f\n". cub010);
			printf("cube(1,1,0) = %f\n". cub110);
			printf("cube(0,0,1) = %f\n". cub001);
			printf("cube(1,0,1) = %f\n". cub101);
			printf("cube(0,1,1) = %f\n". cub011);
			printf("cube(1,1,1) = %f\n". cub111);
#endif

			// Calculating affinity energy
			partialE1 = cub000 * weight000 +
						cub100 * weight100 +
						cub010 * weight010 +
						cub110 * weight110 +
						cub001 * weight001 +
						cub101 * weight101 +
						cub011 * weight011 +
						cub111 * weight111;

#ifdef (PRINT_ALL)
			printf("interpolated energy partialE1 = %f\n\n", partialE1);
#endif

			// -------------------------------------------------------------------
			// TODO: Deltas dx, dy, dz are already normalized
			// (in host/src/getparameters.cpp) in AD-GPU
			// The correspondance between vertices in the xyz axes is:
			// 0, 1, 2, 3, 4, 5, 6, 7 and 000, 100, 010, 001, 101, 110, 011, 111
			// -------------------------------------------------------------------
			/*
				deltas: (x-x0)/(x1-x0), (y-y0...
				vertices: (000, 100, 010, 001, 101, 110, 011, 111)
			  Z
			  '
			  3 - - - - 6
			 /.        /|
			4 - - - - 7 |
			| '       | |
			| 0 - - - + 2 -- Y
			'/        |/
			1 - - - - 5
		       /
		      X
			*/

			// See detailed decomposition in original AD-GPU

			// -------------------------------------------------------------------
			// Calculating gradients (forces) corresponding to
			// "atype" intermolecular energy
			// -------------------------------------------------------------------

			// Vector in x-direction
			gradient_inter_x[atom_id] += omdz * (omdy * (cub100 - cub000) + dy * (cub110 - cub010)) +
										   dz * (omdy * (cub101 - cub001) + dy * (cub111 - cub011));

			// Vector in y-direction
			gradient_inter_y[atom_id] += omdz * (omdx * (cub010 - cub000) + dx * (cub110 - cub100)) +
										   dz * (omdx * (cub011 - cub001) + dx * (cub111 - cub101));

			// Vector in z-direction
			gradient_inter_z[atom_id] += omdy * (omdx * (cub001 - cub000) + dx * (cub101 - cub100)) +
										   dy * (omdx * (cub011 - cub010) + dx * (cub111 - cub110));

			// Energy contribution of the electrostatic grid
			cub000 = (*IE_Fg_2)[iz  ][iy  ][ix  ];
			cub100 = (*IE_Fg_2)[iz  ][iy  ][ix+1];
			cub010 = (*IE_Fg_2)[iz  ][iy+1][ix  ];
			cub110 = (*IE_Fg_2)[iz  ][iy+1][ix+1];
			cub001 = (*IE_Fg_2)[iz+1][iy  ][ix  ];
			cub101 = (*IE_Fg_2)[iz+1][iy  ][ix+1];
			cub011 = (*IE_Fg_2)[iz+1][iy+1][ix  ];
			cub111 = (*IE_Fg_2)[iz+1][iy+1][ix+1];

#ifdef (PRINT_ALL)
			printf("Interpolation of electrostatic map:\n");
			printf("cube(0,0,0) = %f\n". cub000);
			printf("cube(1,0,0) = %f\n". cub100);
			printf("cube(0,1,0) = %f\n". cub010);
			printf("cube(1,1,0) = %f\n". cub110);
			printf("cube(0,0,1) = %f\n". cub001);
			printf("cube(1,0,1) = %f\n". cub101);
			printf("cube(0,1,1) = %f\n". cub011);
			printf("cube(1,1,1) = %f\n". cub111);
#endif

			partialE2 = q * (
						cub000 * weight000 +
						cub100 * weight100 +
						cub010 * weight010 +
						cub110 * weight110 +
						cub001 * weight001 +
						cub101 * weight101 +
						cub011 * weight011 +
						cub111 * weight111);

#ifdef (PRINT_ALL)
			printf("q =%f, interpolated energy partialE2 = %f\n\n", q, partialE2);
#endif

			// -------------------------------------------------------------------
			// Calculating gradients (forces) corresponding to
			// "elec" intermolecular energy
			// -------------------------------------------------------------------

			// Vector in x-direction
			gradient_inter_x[atom_id] += q * (omdz * (omdy * (cub100 - cub000) + dy * (cub110 - cub010)) +
										        dz * (omdy * (cub101 - cub001) + dy * (cub111 - cub011)));

			// Vector in y-direction
			gradient_inter_y[atom_id] += q * (omdz * (omdx * (cub010 - cub000) + dx * (cub110 - cub100)) +
										        dz * (omdx * (cub011 - cub001) + dx * (cub111 - cub101)));

			// Vector in z-direction
			gradient_inter_z[atom_id] += q * (omdy * (omdx * (cub001 - cub000) + dx * (cub101 - cub100)) +
										        dy * (omdx * (cub011 - cub010) + dx * (cub111 - cub110)));

			// Energy contribution of the desolvation grid
			cub000 = (*IE_Fg_3)[iz  ][iy  ][ix  ];
			cub100 = (*IE_Fg_3)[iz  ][iy  ][ix+1];
			cub010 = (*IE_Fg_3)[iz  ][iy+1][ix  ];
			cub110 = (*IE_Fg_3)[iz  ][iy+1][ix+1];
			cub001 = (*IE_Fg_3)[iz+1][iy  ][ix  ];
			cub101 = (*IE_Fg_3)[iz+1][iy  ][ix+1];
			cub011 = (*IE_Fg_3)[iz+1][iy+1][ix  ];
			cub111 = (*IE_Fg_3)[iz+1][iy+1][ix+1];

#ifdef (PRINT_ALL)
			printf("Interpolation of desolvation map:\n");
			printf("cube(0,0,0) = %f\n". cub000);
			printf("cube(1,0,0) = %f\n". cub100);
			printf("cube(0,1,0) = %f\n". cub010);
			printf("cube(1,1,0) = %f\n". cub110);
			printf("cube(0,0,1) = %f\n". cub001);
			printf("cube(1,0,1) = %f\n". cub101);
			printf("cube(0,1,1) = %f\n". cub011);
			printf("cube(1,1,1) = %f\n". cub111);
#endif

			float fabsf_q = fabsf(q);

			partialE3 = fabsf_q * (
						cub000 * weight000 +
						cub100 * weight100 +
						cub010 * weight010 +
						cub110 * weight110 +
						cub001 * weight001 +
						cub101 * weight101 +
						cub011 * weight011 +
						cub111 * weight111);

#ifdef (PRINT_ALL)
			printf("fabsf(q) =%f, interpolated energy partialE3 = %f\n\n", fabsf_q, partialE3);
#endif

			// -------------------------------------------------------------------
			// Calculating gradients (forces) corresponding to
			// "dsol" intermolecular energy
			// -------------------------------------------------------------------

			// Vector in x-direction
			gradient_inter_x[atom_id] += fabsf_q * (omdz * (omdy * (cub100 - cub000) + dy * (cub110 - cub010)) +
													  dz * (omdy * (cub101 - cub001) + dy * (cub111 - cub011)));

			// Vector in y-direction
			gradient_inter_y[atom_id] += fabsf_q * (omdz * (omdx * (cub010 - cub000) + dx * (cub110 - cub100)) +
													  dz * (omdx * (cub011 - cub001) + dx * (cub111 - cub101)));

			// Vector in z-direction
			gradient_inter_z[atom_id] += fabsf_q * (omdy * (omdx * (cub001 - cub000) + dx * (cub101 - cub100)) +
													  dy * (omdx * (cub011 - cub010) + dx * (cub111 - cub110)));

		} // End if

		*final_interE += partialE1 + partialE2 + partialE3; // TODO: eventually will use final_interE[j]

	} // End for (uint atom_id = 0 ...)

	float delta_distance = 0.5f*DockConst_smooth;

	// For each intramolecular atom contributor pair
	for (uint contributor_counter = 0; contributor_counter < DockConst_num_of_intraE_contributors; contributor_counter++)
	{
		// The gradient contribution of each contributing atomic pair
		float priv_gradient_per_intracontributor = 0.0f;

		int atom1_id = IA_intraE_contributors[3*contributor_counter];
		int atom2_id = IA_intraE_contributors[3*contributor_counter + 1];
		int is_H_bond = IA_intraE_contributors[3*contributor_counter + 2];
		uint hbond = (is_H_bond == 1)? 1:0; // evaluates to 1 in case of H-bond, 0 otherwise	// TODO: apply to Solis-Wets

		// Getting types ids
		int atom1_typeid = IA_IE_atom_types[atom1_id];
		int atom2_typeid = IA_IE_atom_types[atom2_id];

		int atom1_type_vdw_hb = IA_atom1_types_reqm[atom1_typeid];
		int atom2_type_vdw_hb = IA_atom2_types_reqm[atom2_typeid];

		// Calculation of delta distance moved outside this loop
		// TODO: test the same for Solis-Wets

		// Getting optimum pair distance (opt_distance) from reqm and reqm_hbond
		// reqm: equilibrium internuclear separation
		//       (sum of the vdW radii of two like atoms (A)) in the case of vdW
		// reqm_hbond: equilibrium internuclear separation
		// 	 (sum of the vdW radii of two like atoms (A)) in the case of hbond
		float opt_distance;

		if (hbond) {	// H-bond
			opt_distance = IA_reqm_hbond[atom1_type_vdw_hb] + IA_reqm_hbond[atom2_type_vdw_hb];
		} else {	// Van der Waals
			opt_distance = 0.5f*(IA_reqm[atom1_type_vdw_hb] + IA_reqm[atom2_type_vdw_hb]);
		}

		float vdW_const1 = IA_VWpars_AC[atom1_typeid * DockConst_num_of_atypes + atom2_typeid];
		float vdW_const2 = IA_VWpars_BD[atom1_typeid * DockConst_num_of_atypes + atom2_typeid];

		float desolv_const = ((IA_dspars_S[atom1_typeid] + DockConst_qasp * esa_fabs(IA_IE_atom_charges[atom1_id])) * IA_dspars_V[atom2_typeid] +
				      		  (IA_dspars_S[atom2_typeid] + DockConst_qasp * esa_fabs(IA_IE_atom_charges[atom2_id])) * IA_dspars_V[atom1_typeid]
							 ) * DockConst_coeff_desolv;
		float elec_const = DockConst_coeff_elec * IA_IE_atom_charges[atom1_id] * IA_IE_atom_charges[atom2_id];

		// TODO: add support for flexible rings

		// TODO: eventually add loop for vectorization
		//for (uint j = 0; j < DockConst_pop_size; j++) {
			float subx = local_coords_x[atom1_id] - local_coords_x[atom2_id];
			float suby = local_coords_y[atom1_id] - local_coords_y[atom2_id];
			float subz = local_coords_z[atom1_id] - local_coords_z[atom2_id];

			float atomic_distance = esa_sqrt(subx*subx + suby*suby + subz*subz)*DockConst_grid_spacing;

#ifdef (PRINT_ALL)
			printf("\nContrib %u: atoms %u and %u, distance: %f\n", contributor_counter, atom1_id+1, atom2_id+1, atomic_distance);
#endif

			float partialIAE1 = 0.0f;
			float partialIAE2 = 0.0f;
			float partialIAE3 = 0.0f;
			float partialIAE4 = 0.0f;

			// Getting smoothed_distance = function(atomic_distance, opt_distance)
			float smoothed_distance;

			if (atomic_distance <= (opt_distance - delta_distance)) {
				smoothed_distance = atomic_distance + delta_distance;
			}
			else if (atomic_distance < (opt_distance + delta_distance)) {
				smoothed_distance = opt_distance;
			}
			else { // else if (atomic_distance >= (opt_distance + delta_distance))
				smoothed_distance = atomic_distance - delta_distance;
			}

			float distance_pow_2  = atomic_distance * atomic_distance;

			float inverse_smoothed_distance = 1.0f / smoothed_distance;
			float inverse_smoothed_distance_pow_2  = inverse_smoothed_distance * inverse_smoothed_distance;
			float inverse_smoothed_distance_pow_4  = inverse_smoothed_distance_pow_2 * inverse_smoothed_distance_pow_2;
			float inverse_smoothed_distance_pow_6  = inverse_smoothed_distance_pow_4 * inverse_smoothed_distance_pow_2;
			float inverse_smoothed_distance_pow_12 = inverse_smoothed_distance_pow_6 * inverse_smoothed_distance_pow_6;	// TODO: fix it in Solis-wets

			// Calculating energy contributions
			// Cuttoff1: internuclear-distance at 8A only for vdw and hbond.
			if (atomic_distance < 8.0f) {

				partialIAE1 = vdW_const1 * inverse_smoothed_distance_pow_12; // TODO: do the same for Solis-Wets
				float partialIAE1_times_12 = 12.0f * partialIAE1;

				float gradient_numerator;

				// Calculating van der Waals / hydrogen bond term
				if (hbond) {	// H-bond
					float inverse_smoothed_distance_pow_10 = inverse_smoothed_distance_pow_6 * inverse_smoothed_distance_pow_4;
					partialIAE2 = vdW_const2 * inverse_smoothed_distance_pow_10;
					gradient_numerator = 10.0f * partialIAE2 - partialIAE1_times_12;
				}
				else {	// Van der Waals
					partialIAE2 = vdW_const2 * inverse_smoothed_distance_pow_6;
					gradient_numerator =  6.0f * partialIAE2 - partialIAE1_times_12;
				}

				priv_gradient_per_intracontributor += gradient_numerator * inverse_smoothed_distance;

			} // End if cuttoff1 - internuclear-distance at 8A

			// Calculating energy contributions
			// Cuttoff2: internuclear-distance at 20.48A only for el and sol.
			if (atomic_distance < 20.48f) {
				float term_partialE3 = atomic_distance *
					(DIEL_A + (DIEL_B / (1.0f + DIEL_K * esa_expf0(-DIEL_B_TIMES_H * atomic_distance))));
				float term_inv_partialE3 = (1.0f / term_partialE3);

				// Calculating electrostatic term
				partialIAE3 = elec_const * term_inv_partialE3;

				// Calculating desolvation term
				partialIAE4 = desolv_const * esa_expf0(-0.03858025f * distance_pow_2);

				//priv_gradient_per_intracontributor += ;

			} // End if cuttoff2 - internuclear-distance at 20.48A

		//} // End for (uint j = 0 ...)

		*final_intraE += partialIAE1 - partialIAE2 + partialIAE3 + partialIAE4;

	} // End for (uint contributor_counter = 0 ...)

}
