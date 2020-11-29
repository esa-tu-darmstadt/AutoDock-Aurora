#include "auxiliary.c"

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
void perform_ls(
			ushort				DockConst_max_num_of_iters,
			float				DockConst_rho_lower_bound,
			float				DockConst_base_dmov_mul_sqrt3,
			uchar				DockConst_num_of_genes,
			float				DockConst_base_dang_mul_sqrt3,
			uchar				DockConst_cons_limit,

			float*	restrict	in_out_genotype,
			float*	restrict	in_out_energy,
			uint*	restrict	out_eval,
			uint*				dockpars_prng_states,
	// pc
	const	int*	restrict	PC_rotlist,
	const	float*	restrict	PC_ref_coords_x,
	const	float*	restrict	PC_ref_coords_y,
	const	float*	restrict	PC_ref_coords_z,
	const	float*	restrict	PC_rotbonds_moving_vectors,
	const	float*	restrict	PC_rotbonds_unit_vectors,
	const	float*	restrict	PC_ref_orientation_quats,
			uint				DockConst_rotbondlist_length,
			uint				Host_RunId,
	// ia
	const	float*	restrict	IA_IE_atom_charges,
	const	int*	restrict	IA_IE_atom_types,
	const	int*	restrict	IA_intraE_contributors,
	const	float*	restrict	IA_reqm,
	const	float*	restrict	IA_reqm_hbond,
	const	uint*	restrict	IA_atom1_types_reqm,
	const	uint*	restrict	IA_atom2_types_reqm,
	const	float*	restrict	IA_VWpars_AC,
	const	float*	restrict	IA_VWpars_BD,
	const	float*	restrict	IA_dspars_S,
	const	float*	restrict	IA_dspars_V,
			float				DockConst_smooth,
			uint				DockConst_num_of_intraE_contributors,
			float				DockConst_grid_spacing,
			uint				DockConst_num_of_atypes,
			float				DockConst_coeff_elec,
			float				DockConst_qasp,
			float				DockConst_coeff_desolv,
	// ie
	const	float*	restrict	IE_Fgrids,
			uchar				DockConst_g1,
			uint				DockConst_g2,
			uint				DockConst_g3,
			uchar				DockConst_num_of_atoms,
			float				DockConst_gridsize_x_minus1,
			float				DockConst_gridsize_y_minus1,
			float				DockConst_gridsize_z_minus1,
			uint				Host_mul_tmp2,
			uint				Host_mul_tmp3
)
{	
	#if defined (PRINT_ALL_LS) 
	printf("\n");
	printf("Starting <local search> ... \n");
	printf("\n");
	printf("LS: DockConst_max_num_of_iters: %u\n",		DockConst_max_num_of_iters);
	printf("LS: DockConst_rho_lower_bound: %f\n",      	DockConst_rho_lower_bound);
	printf("LS: DockConst_base_dmov_mul_sqrt3: %f\n",  	DockConst_base_dmov_mul_sqrt3);
	printf("LS: DockConst_num_of_genes: %u\n",  	   	DockConst_num_of_genes);
	printf("LS: DockConst_base_dang_mul_sqrt3: %f\n",  	DockConst_base_dang_mul_sqrt3);
	printf("LS: DockConst_cons_limit: %u\n",           	DockConst_cons_limit);
	#endif

	float local_coords_x[MAX_NUM_OF_ATOMS];
	float local_coords_y[MAX_NUM_OF_ATOMS];
	float local_coords_z[MAX_NUM_OF_ATOMS];

	// Reading incoming genotype and energy
	float current_energy = *in_out_energy;
	float genotype[ACTUAL_GENOTYPE_LENGTH];
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
		genotype[i] = in_out_genotype[i];
	}
	
	float rho = 1.0f;
	ushort iteration_cnt = 0;
	uchar  cons_succ     = 0;
	uchar  cons_fail     = 0;
	uint   LS_eval       = 0;
	boolean positive_direction = True;

	// Performing local search
	while ((iteration_cnt < DockConst_max_num_of_iters) && (rho > DockConst_rho_lower_bound)) {
		// -----------------------------------------------
		// Exit condition is groups here. It allows pipelining
		if (positive_direction == True) {
			if (cons_succ >= DockConst_cons_limit) {
				rho = LS_EXP_FACTOR * rho;
				cons_fail = 0;
				cons_succ = 0;
			}
			else if (cons_fail >= DockConst_cons_limit) {
				rho = LS_CONT_FACTOR * rho;
				cons_fail = 0;
				cons_succ = 0;
			}
			iteration_cnt++;
		}

		#if defined (PRINT_ALL)
		printf("LS positive?: %u, iteration_cnt: %u, rho: %f, limit rho: %f\n",
		positive_direction, iteration_cnt, rho, DockConst_rho_lower_bound);
		#endif
		// -----------------------------------------------

		float entity_possible_new_genotype[ACTUAL_GENOTYPE_LENGTH];
		float genotype_bias[ACTUAL_GENOTYPE_LENGTH];
		float deviate_plus_bias[ACTUAL_GENOTYPE_LENGTH];
		float deviate_minus_bias[ACTUAL_GENOTYPE_LENGTH];

		// Generating new random deviate
		// rho is the deviation of the uniform distribution
		for (uint i = 0; i < DockConst_num_of_genes; i++) {
			// TODO: FIX INDEXES
			float tmp_prng = randf(dockpars_prng_states);

			// tmp1 is genotype_deviate
			float tmp1 = rho * (2.0f * tmp_prng - 1.0f);

			if (i < 3) {
				tmp1 = tmp1 * DockConst_base_dmov_mul_sqrt3;
			}
			else {
				tmp1 = tmp1 * DockConst_base_dang_mul_sqrt3;
			}

			float deviate = 0.4f * tmp1;

			// tmp2 is the addition: genotype_deviate + genotype_bias
			// tmp3 is entity_possible_new_genotype
			float tmp_bias = (iteration_cnt == 1)? 0.0f : genotype_bias[i];
			float bias = 0.6f * tmp_bias;

			deviate_plus_bias[i] = deviate + bias;
			deviate_minus_bias[i] = deviate - bias;

			float tmp2 = tmp1 + tmp_bias;
			float tmp3 = (positive_direction == True)? (genotype [i] + tmp2): (genotype [i] - tmp2);

			if (i > 2) {
				if (i == 4) {
					tmp3 = map_angle_180(tmp3);
				}
				else {
					tmp3 = map_angle_360(tmp3);
				}
			}

			entity_possible_new_genotype[i] = tmp3;
		}

		// TODO: CALC ENERGY
		float energy_ia_ls;
		float energy_ie_ls;
		calc_pc(
			PC_rotlist,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			DockConst_rotbondlist_length,
			DockConst_num_of_genes,
			Host_RunId,
			entity_possible_new_genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z
		);
		energy_ia(
			IA_IE_atom_charges,
			IA_IE_atom_types,
			IA_intraE_contributors,
			IA_reqm,
			IA_reqm_hbond,
			IA_atom1_types_reqm,
			IA_atom2_types_reqm,
			IA_VWpars_AC,
			IA_VWpars_BD,
			IA_dspars_S,
			IA_dspars_V,
			DockConst_smooth,
			DockConst_num_of_intraE_contributors,
			DockConst_grid_spacing,
			DockConst_num_of_atypes,
			DockConst_coeff_elec,
			DockConst_qasp,
			DockConst_coeff_desolv,
			&energy_ia_ls,
			local_coords_x,
			local_coords_y,
			local_coords_z
		);
		energy_ie(
			IE_Fgrids,
			IA_IE_atom_charges,
			IA_IE_atom_types,
			DockConst_g1,
			DockConst_g2,
			DockConst_g3,
			DockConst_num_of_atoms,
			DockConst_gridsize_x_minus1,
			DockConst_gridsize_y_minus1,
			DockConst_gridsize_z_minus1,
			Host_mul_tmp2,
			Host_mul_tmp3,
			&energy_ie_ls,
			local_coords_x,
			local_coords_y,
			local_coords_z
		);
		float candidate_energy = energy_ia_ls + energy_ie_ls;

		// Updating LS energy-evaluation count
		LS_eval++;

		if (candidate_energy < current_energy) {
			// Updating offspring_genotype & genotype_bias

			for (uint i = 0; i < DockConst_num_of_genes; i++) {
				genotype_bias[i] = (positive_direction == True) ? deviate_plus_bias[i] : deviate_minus_bias[i];
				genotype[i] = entity_possible_new_genotype[i];
			}

			current_energy = candidate_energy;
			cons_succ++;
			cons_fail = 0;
			positive_direction = True;
		}
		else {
			// Updating (halving) genotype_bias

			for (uint i = 0; i < DockConst_num_of_genes; i++) {
				genotype_bias[i] = (iteration_cnt == 1)? 0.0f: (0.5f*genotype_bias[i]);
			}

			if (positive_direction == False) {
				cons_fail++;
				cons_succ = 0;
			}
			positive_direction = (positive_direction == True) ? False: True;
		}

	} // end of while (iteration_cnt) && (rho)

	// Writing resulting number of energy evals performed in LS
	*out_eval = LS_eval;

	// Writing resulting genotype and energy
	*in_out_energy = current_energy;
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
		in_out_genotype[i] = genotype[i];
	}
	
	#if defined (PRINT_ALL_IE) 
	printf("\n");
	printf("Finishing <local search>\n");
	printf("\n");
	#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
