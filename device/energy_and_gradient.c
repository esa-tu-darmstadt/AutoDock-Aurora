#include "auxiliary.h"

void energy_and_gradient (
			float*				genotype,
	const	uchar               DockConst_num_of_genes, // ADGPU defines it as int
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
			float*				gradient_genotype,
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

	// Initializing
	// TODO: make sure this is strictly necessary
	// (caller may be doing the same, but maybe not redundant)
	// TODO: merge with for-loop below
	for (uint atom_id = 0; atom_id < DockConst_num_of_atoms; atom_id++) {
		gradient_inter_x[atom_id] = 0.0f;
		gradient_inter_y[atom_id] = 0.0f;
		gradient_inter_z[atom_id] = 0.0f;

		gradient_intra_x[atom_id] = 0.0f;
		gradient_intra_y[atom_id] = 0.0f;
		gradient_intra_z[atom_id] = 0.0f;
	}

	for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
		gradient_genotype[gene_cnt] = 0.0f;
	}

	// ================================================
	// CALCULATING INTERMOLECULAR ENERGY & GRADIENTS
	// ================================================
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
			float cub000, cub001, cub010, cub011;
			float cub100, cub101, cub110, cub111;

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
			partialE1 = cub000 * weight000 + cub100 * weight100 +
						cub010 * weight010 + cub110 * weight110 +
						cub001 * weight001 + cub101 * weight101 +
						cub011 * weight011 + cub111 * weight111;

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
						cub000 * weight000 + cub100 * weight100 +
						cub010 * weight010 + cub110 * weight110 +
						cub001 * weight001 + cub101 * weight101 +
						cub011 * weight011 + cub111 * weight111);

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
						cub000 * weight000 + cub100 * weight100 +
						cub010 * weight010 + cub110 * weight110 +
						cub001 * weight001 + cub101 * weight101 +
						cub011 * weight011 + cub111 * weight111);

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

	// ================================================
	// CALCULATING INTRAMOLECULAR ENERGY & GRADIENTS
	// ================================================

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

			// Calculating atomic distance
			float dist = esa_sqrt(subx*subx + suby*suby + subz*subz);
			float atomic_distance = dist * DockConst_grid_spacing;

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
			float inverse_smoothed_distance_pow_13 = inverse_smoothed_distance_pow_12 * inverse_smoothed_distance;

			// Calculating energy contributions
			// Cuttoff1: internuclear-distance at 8A only for vdw and hbond.
			if (atomic_distance < 8.0f) {

				partialIAE1 = vdW_const1 * inverse_smoothed_distance_pow_12; // TODO: do the same for Solis-Wets
				float partialIAE1_times_12 = 12.0f * partialIAE1;

				priv_gradient_per_intracontributor += (-12.0f * vdW_const1) * inverse_smoothed_distance_pow_13;

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
				float tmp_exp0 = esa_expf0(DIEL_B_TIMES_H * atomic_distance);
				float inv_tmp_exp0 = 1.0f / tmp_exp0;

				float term_partialE3 = atomic_distance * (DIEL_A + (DIEL_B / (1.0f + DIEL_K * inv_tmp_exp0)));
				float term_inv_partialE3 = (1.0f / term_partialE3);

				// Calculating electrostatic term
				partialIAE3 = elec_const * term_inv_partialE3;

				// http://www.wolframalpha.com/input/?i=1%2F(x*(A%2B(B%2F(1%2BK*exp(-h*B*x)))))
				float tmp_exp1 = tmp_exp0 + DIEL_K;
				float tmp_exp2 = tmp_exp0 * (DIEL_B_TIMES_H_TIMES_K * atomic_distance + tmp_exp0 + DIEL_K);
				float upper = DIEL_A * tmp_exp1 * tmp_exp1 + DIEL_B * tmp_exp2;

				float tmp_exp3 = DIEL_A * (tmp_exp0 + DIEL_K) + DIEL_B * tmp_exp0;
				float lower = distance_pow_2 * tmp_exp3 * tmp_exp3;

				priv_gradient_per_intracontributor += -elec_const * (upper/lower);

				// Calculating desolvation term
				float tmp_exp4 = esa_expf0(-0.03858025f * distance_pow_2);
				partialIAE4 = desolv_const * tmp_exp4;

				priv_gradient_per_intracontributor += -0.077160f * atomic_distance * partialIAE4;

			} // End if cuttoff2 - internuclear-distance at 20.48A

		//} // End for (uint j = 0 ...)

		*final_intraE += partialIAE1 - partialIAE2 + partialIAE3 + partialIAE4;

		// Decomposing "priv_gradient_per_intracontributor"
		// into the contribution of each atom of the pair.
		// Distances in Ansgtroms of vector that goes from "atom1_id"-to-"atom2_id".
		// Therefore, -subx, -suby. -subz are used.
		float subx_div_dist = -subx / dist;
		float suby_div_dist = -suby / dist;
		float subz_div_dist = -subz / dist;

		float priv_intra_gradient_x = priv_gradient_per_intracontributor * subx_div_dist;
		float priv_intra_gradient_y = priv_gradient_per_intracontributor * suby_div_dist;
		float priv_intra_gradient_z = priv_gradient_per_intracontributor * subz_div_dist;

		// Calculating gradients in xyz components
		// Gradients for both atoms in a single contributor pair
		// have the same magnitude, but opposite directions

		// TODO: confirm no need for atomic ops (because this is not SIMT as in AD-GPU)
		gradient_intra_x[atom1_id] = gradient_intra_x[atom1_id] - priv_intra_gradient_x;
		gradient_intra_y[atom1_id] = gradient_intra_y[atom1_id] - priv_intra_gradient_y;
		gradient_intra_z[atom1_id] = gradient_intra_z[atom1_id] - priv_intra_gradient_z;

		gradient_intra_x[atom1_id] = gradient_intra_x[atom2_id] + priv_intra_gradient_x;
		gradient_intra_y[atom1_id] = gradient_intra_y[atom2_id] + priv_intra_gradient_y;
		gradient_intra_z[atom1_id] = gradient_intra_z[atom2_id] + priv_intra_gradient_z;
	} // End for (uint contributor_counter = 0 ...)

	// ================================================
	// ACCUMMULATING INTER & INTRAMOLECULAR GRADIENTS
	// ================================================
	for (uint atom_cnt = 0; atom_cnt < DockConst_num_of_atoms; atom_cnt++) {

		// Grids were calculated in the grid space,
		// so they have to be put back in Angstrom

		// Intramolecular gradients were already in Angstrom,
		// no scaling for them is required

		 // TODO: compute ind in host and then pass it
		gradient_inter_x[atom_id] = gradient_inter_x[atom_id] / DockConst_grid_spacing;
		gradient_inter_y[atom_id] = gradient_inter_y[atom_id] / DockConst_grid_spacing;
		gradient_inter_z[atom_id] = gradient_inter_z[atom_id] / DockConst_grid_spacing;

		// Reusing "gradient_inter_*" for total gradient (inter + intra)
		gradient_inter_x[atom_id] += gradient_intra_x[atom_id];
		gradient_inter_y[atom_id] += gradient_intra_y[atom_id];
		gradient_inter_z[atom_id] += gradient_intra_z[atom_id];
	}

	// ================================================
	// OBTAINING TRANSLATION-RELATED GRADIENTS
	// ================================================
	for (uint atom_id = 0; atom_id < DockConst_num_of_atoms; atom_id++) {
		// Accummulating "gradient_inter_*"
		gradient_genotype[0] += gradient_inter_x[atom_id]; // gradient for gene 0: gene x
		gradient_genotype[1] += gradient_inter_y[atom_id]; // gradient for gene 1: gene y
		gradient_genotype[2] += gradient_inter_z[atom_id]; // gradient for gene 2: gene z
	}

	// Scaling gradient for translational genes as
	// their corresponding gradients were calculated in the space
	// where these genes are in Angstrom, but AD-GPU translational genes are in grids
	gradient_genotype[0] *= DockConst_grid_spacing;
	gradient_genotype[1] *= DockConst_grid_spacing;
	gradient_genotype[2] *= DockConst_grid_spacing;

	// ================================================
	// OBTAINING ROTATION-RELATED GRADIENTS

	// Transform gradient_inter_{x|y|z}
	// into local_gradients[i] (with four quaternion genes)

	// Transform local_gradients[i] (with four quaternion genes)
	// into local_gradients[i] (with three genes)
	// ================================================

	float torque_rot_x = 0.0f;
	float torque_rot_y = 0.0f;
	float torque_rot_z = 0.0f;

	// Variable holding the center of rotation
	// In getparameters.cpp, it indicates translation genes are
	// in grid spacing (instead of Angstrom)
	float about_x = genotype[0];
	float about_y = genotype[1];
	float about_z = genotype[2];

	// Temporal variable for calculating translation differences.
	// They are converted back to Angstrom here
	float r_x, r_y; r_z;

	for (uint atom_id = 0; atom_id < DockConst_num_of_atoms; atom_id++) {
		r_x = (local_coords_x[atom_id] - about_x) * DockConst_grid_spacing;
		r_y = (local_coords_y[atom_id] - about_y) * DockConst_grid_spacing;
		r_z = (local_coords_z[atom_id] - about_z) * DockConst_grid_spacing;

		// Reusing "gradient_inter_*" for total gradient (inter + intra)
		float force_x = gradient_inter_x[atom_id];
		float force_y = gradient_inter_y[atom_id];
		float force_z = gradient_inter_z[atom_id];

		float tmp_x, tmp_y, tmp_z;
		esa_cross3_e_(r_x, r_y, r_z, force_x, force_y, force_z, &tmp_x, &tmp_y, &tmp_z);
		torque_rot_x += tmp_x;
		torque_rot_y += tmp_y;
		torque_rot_z += tmp_z;
	}

	// genes[3:7] = rotation.axisangle_to_q(torque, rad)
	float torque_length = esa_length3_e(torque_rot_x, torque_rot_y, torque_rot_z);

	// Finding quaternion performing
	// infinitesimal rotation around torque axis
	float quat_torque_x, quat_torque_y, quat_torque_z, quat_torque_w;

	float tmp_normal_x, tmp_normal_y, tmp_normal_z;
	esa_normalize3_e_(torque_rot_x, torque_rot_y, torque_rot_z, tmp_normal_x, tmp_normal_y, tmp_normal_z);

	quat_torque_w = COS_HALF_INFINITESIMAL_RADIAN;
	quat_torque_x = tmp_normal_x * SIN_HALF_INFINITESIMAL_RADIAN;
	quat_torque_y = tmp_normal_y * SIN_HALF_INFINITESIMAL_RADIAN;
	quat_torque_z = tmp_normal_z * SIN_HALF_INFINITESIMAL_RADIAN;

	// Converting quaternion gradients into orientation gradients

	// This is where we are in the orientation axis-angle space
	// TODO: check very initial input orientation genes
	float current_phi, current_theta, current_rotangle;
	current_phi 	 = genotype[3];	// phi		(in sexagesimal (DEG) unbounded)
	current_theta 	 = genotype[4];	// theta	(in sexagesimal (DEG) unbounded)
	current_rotangle = genotype[5];	// rotangle	(in sexagesimal (DEG) unbounded)

	current_phi   = map_angle_360 (current_phi);	// phi (in DEG bounded [0, 360])
	current_theta = map_angle_180 (current_theta);	// theta (in DEG bounded [0, 180]) // TODO: should it bounded to [0, 180] ??
	current_rotangle = map_angle_360 (current_rotangle);	// rotangle (in DEG bounded [0, 360])

	current_phi = current_phi * DEG_TO_RAD;	// phi (in GRAD)
	current_theta = current_theta * DEG_TO_RAD;	// theta (in GRAD)
	current_rotangle = current_rotangle * DEG_TO_RAD;	// rotangle (in GRAD)

	int is_theta_gt_pi = (current_theta > PI_FLOAT) ? 1 : 0;

	// This is where we are in the quaternion space
	float current_q_w, current_q_x, current_q_y, current_q_z;

	// Axis of rotation
	float rotaxis_x = sinf(current_theta) * cosf(current_phi);
	float rotaxis_y = sinf(current_theta) * sinf(current_phi);
	float rotaxis_z = cosf(current_theta);

	float ang = current_rotangle * 0.5f;
	current_q_w = cosf(ang);
	current_q_x = rotaxis_x * sinf(ang);
	current_q_y = rotaxis_y * sinf(ang);
	current_q_z = rotaxis_z * sinf(ang);

	// This is where we want to be in the quaternion space
	float target_q_w, target_q_x, target_q_y, target_q_z;

	// target_q = rotation.q_mult(q, current_q)
	// In our terms it means: q_mult(quat_{w|x|y|z}, currrent_q{w|x|y|z})
	//target_q_w = quat_torque_w * current_q_w - quat_torque_x * current_q_x - quat_torque_y * current_q_y - quat_torque_z * current_q_z;
	//target_q_x = quat_torque_w * current_q_x + quat_torque_x * current_q_w + quat_torque_y * current_q_z - quat_torque_z * current_q_y;
	//target_q_y = quat_torque_w * current_q_y + quat_torque_y * current_q_w + quat_torque_z * current_q_x - quat_torque_x * current_q_z;
	//target_q_z = quat_torque_w * current_q_z + quat_torque_z * current_q_w + quat_torque_x * current_q_y - quat_torque_y * current_q_x;

	target_q_w = esa_dot4_e(quat_torque_w, quat_torque_x, quat_torque_y, quat_torque_z,
							current_q_w, -current_q_x, -current_q_y, -current_q_z);
	target_q_x = esa_dot4_e(quat_torque_w, quat_torque_x, quat_torque_y, quat_torque_z,
							current_q_x, current_q_w, current_q_z, -current_q_y);
	target_q_y = esa_dot4_e(quat_torque_w, quat_torque_x, quat_torque_y, quat_torque_z,
							current_q_y, -current_q_z, current_q_w, current_q_x);
	target_q_z = esa_dot4_e(quat_torque_w, quat_torque_x, quat_torque_y, quat_torque_z,
							current_q_z, current_q_y, -current_q_x, current_q_w);

	// This is where we want to be in the orientation axis-angle space
	float target_phi, target_theta, target_rotangle;

	// target_oclacube = quaternion_to_oclacube(target_q, theta_gt_pi)
	// In our terms it means: quaternion_to_oclacube(target_q{w|x|y|z}, theta_gt_pi)
	ang = acosf(target_q_w);	// TODO: make sure single precision function works!
	target_rotangle = 2.0f * ang;

	float inv_sin_ang = 1.0f / (sinf(ang));
	rotaxis_x = target_q_x * inv_sin_ang;
	rotaxis_y = target_q_y * inv_sin_ang;
	rotaxis_z = target_q_z * inv_sin_ang;

	target_theta = acosf(rotaxis_z); // TODO: make sure single precision function works!

	if (is_theta_gt_pi == 0) { // false
		target_phi = fmodf((atan2f(rotaxis_y, rotaxis_x) + PI_TIMES_2), PI_TIMES_2); // TODO: check if fmod is supported, https://sleef.org/purec.xhtml
	} else {
		target_phi = fmodf((atan2f(-rotaxis_y, -rotaxis_x) + PI_TIMES_2), PI_TIMES_2); // TODO: check if fmod is supported, https://sleef.org/purec.xhtml
		target_theta = PI_TIMES_2 - target_theta;
	}

	// The infinitesimal rotation produces an infinitesimal displacement.
	// This is to guarantee that the direction of the displacement is not distorted.
	float orientation_scaling = torque_length * INV_INFINITESIMAL_RADIAN;

	// Derivatives (or gradients)
	float grad_phi = orientation_scaling * (fmodf(target_phi - current_phi + PI_FLOAT, PI_TIMES_2) - PI_FLOAT);
	float grad_theta = orientation_scaling * (fmodf(target_theta - current_theta + PI_FLOAT, PI_TIMES_2) - PI_FLOAT);
	float grad_rotangle = orientation_scaling * (fmodf(target_rotangle - current_rotangle + PI_FLOAT, PI_TIMES_2) - PI_FLOAT);

	// Derivatives corrections.
	// Constant array have 1000 elements.
	// Each array spans approximatedly from 0.0 to 2*PI.
	// The distance between each x-point (angle_delta) is 2*PI/1000.
	const float angle_delta = 0.00628353f;
	const float inv_angle_delta = 159.154943;

	// Correcting theta gradients interpolating
	// values from correction lookup tables
	// (X0, Y0) and (X1, Y1) are known points.
	// How to find the Y value in the straight line between Y0 and Y1,
	// corresponding to a certain X?
	/*
			| dependence_on_theta_const
			| dependence_on_rotangle_const
			|
			|
			|                        Y1
			|
			|             Y=?
			|    Y0
			|_________________________________ angle_const
			     X0  		X	 	 X1
	*/

	// Finding the index position of "grad_delta" in the "angle_const" array
	uint index_theta = floorf((current_theta - GRAD_angle[0]) * inv_angle_delta);
	uint index_rotangle = floorf((current_rotangle - GRAD_angle[0]) * inv_angle_delta);

	// Interpolating theta values
	// X0 -> index - 1
	// X1 -> index + 1
	// Expressed as weighted average
	// Y = [Y0 * ((X1 - X) / (X1-X0))] +  [Y1 * ((X - X0) / (X1-X0))]
	// Simplified (fewer terms)
	// Y = [Y0 * (X1 - X) + Y1 * (X - X0)] / (X1 - X0)
	// Taking advantage of constants
	// Y = [Y0 * (X1 - X) + Y1 * (X - X0)] * inv_angle_delta
	float X0_theta, Y0_theta;
	float X1_theta, Y1_theta;
	float X_theta;
	float dependance_on_theta;	// Y = dependance_on_theta
	X_theta = current_theta;

	// Using interpolation on out-of-bounds elements results in hang
	if (index_theta <= 0) {
		dependance_on_theta = GRAD_dependence_on_theta[0];

	} else if (index_theta >= 999) {
		dependance_on_theta = GRAD_dependence_on_theta[999];
	} else {
		X0_theta = GRAD_angle[index_theta];
		Y0_theta = GRAD_dependence_on_theta[index_theta];
		X1_theta = GRAD_angle[index_theta + 1];
		Y1_theta = GRAD_dependence_on_theta[index_theta + 1];
		// TODO: make this statement is in the correct place
		dependance_on_theta = (Y0_theta * (X1_theta - X_theta) + Y1_theta * (X_theta - X0_theta)) * inv_angle_delta;
	}

	// Interpolating rotangle values
	float X0_rotangle, Y0_rotangle;
	float X1_rotangle, Y1_rotangle;
	float X_rotangle;
	float dependance_on_rotangle; // Y = dependance_on_rotangle
	X_rotangle = current_rotangle;

	// Using interpolation on previous and/or next elements results in hang
	if (index_rotangle <= 0) {
		dependance_on_rotangle = GRAD_dependence_on_rotangle[0];

	} else if (index_rotangle >= 999) {
		dependance_on_rotangle = GRAD_dependence_on_rotangle[999];
	} else {
		X0_rotangle = GRAD_angle[index_rotangle];
		Y0_rotangle = GRAD_dependence_on_rotangle[index_rotangle];
		X1_rotangle = GRAD_angle[index_rotangle + 1];
		Y1_rotangle = GRAD_dependence_on_rotangle[index_rotangle + 1];
		// TODO: make this statement is in the correct place
		dependance_on_rotangle = (Y0_rotangle * (X1_rotangle - X_rotangle) + Y1_rotangle * (X_rotangle - X0_rotangle)) * inv_angle_delta;
	}

	// Setting gradient rotation-related genotypes in cube
	// Multiplying by DEG_TO_RAD to make it uniform to  DEG (see torsion gradients)
	gradient_genotype[3] = (grad_phi / (dependance_on_theta * dependance_on_rotangle)) * DEG_TO_RAD;
	gradient_genotype[4] = (grad_theta / dependance_on_rotangle) * DEG_TO_RAD;
	gradient_genotype[5] = grad_rotangle * DEG_TO_RAD;

	// ================================================
	// OBTAINING TORSION-RELATED GRADIENTS
	// ================================================
	for (uint rotbond_id = 0; rotbond_id < DockConst_num_of_genes - 6; rotbond_id++) {
		// Querying ids of atoms belonging to the rotatable bond in question
		int atom1_id = GRAD_rotbonds[2 * rotbond_id];
		int atom2_id = GRAD_rotbonds[2 * rotbond_id + 1];

		float atomRef_coords_x = local_coords_x[atom1_id];
		float atomRef_coords_y = local_coords_y[atom1_id];
		float atomRef_coords_z = local_coords_z[atom1_id];

		float rotation_unitvec_x = local_coords_x[atom2_id] - local_coords_x[atom1_id];
		float rotation_unitvec_y = local_coords_y[atom2_id] - local_coords_y[atom1_id];
		float rotation_unitvec_z = local_coords_z[atom2_id] - local_coords_z[atom1_id];

		// Torque of torsions
		float torque_tor_x = 0.0f;
		float torque_tor_y = 0.0f;
		float torque_tor_z = 0.0f;

	}
}
