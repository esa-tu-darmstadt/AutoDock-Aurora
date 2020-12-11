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

                        uint			pop_size,
			float*	restrict	in_out_genotype_,
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
	float (*in_out_genotype)[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE] = (float (*)[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE])in_out_genotype_;

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

	float local_coords_x[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_y[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_z[MAX_NUM_OF_ATOMS][MAX_POPSIZE];

	for (uint i = 0; i < DockConst_num_of_atoms; i++) {
		for (uint j = 0; j < pop_size; j++) {
			local_coords_x[i][j] = 0.0f;
			local_coords_y[i][j] = 0.0f;
			local_coords_z[i][j] = 0.0f;
		}
	}

	// Reading incoming genotype and energy
	float current_energy[MAX_POPCOUNT];
	float genotype[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

	for (uint i = 0; i < DockConst_num_of_genes; i++) {
		for (uint j = 0; j < pop_size; j++) {
			genotype[i][j] = in_out_genotype[i][j];
		}
	}
	for (uint j = 0; j < pop_size; j++) {
		current_energy[j] = in_out_energy[j];
	}

        float rho[MAX_POPSIZE];
        uint iteration_cnt[MAX_POPSIZE];
        uint cons_succ[MAX_POPSIZE];
        uint cons_fail[MAX_POPSIZE];
        uint LS_eval[MAX_POPSIZE];
        int positive_direction[MAX_POPSIZE];  // converted from boolean to int
	int ls_is_active[MAX_POPSIZE];  // filter for individuals that are still being iterated
	uint num_active_ls = pop_size;

	int active_pop_size; // counts compressed list of genomes
	uint active_idx[MAX_POPSIZE]; // j index that corresponds to slot in compressed list
	uint active_compr_idx[MAX_POPSIZE]; // compressed index corresponding to active j index
	
	for (uint j = 0; j < pop_size; j++) {
		rho[j]                = 1.0f;
		iteration_cnt[j]      = 0;
		cons_succ[j]          = 0;
		cons_fail[j]          = 0;
		LS_eval[j]            = 0;
		positive_direction[j] = 1;
		ls_is_active[j]       = 1;
	}

	// Performing local search
	//while ((iteration_cnt < DockConst_max_num_of_iters) && (rho > DockConst_rho_lower_bound)) {
	while (num_active_ls > 0) {

		// compressed list of active indices
		active_pop_size = 0;
		for (uint j = 0; j < pop_size; j++) {
			if (ls_is_active[j]) {
				active_idx[active_pop_size] = j;
				active_compr_idx[j] = active_pop_size;
				active_pop_size++;
			}
		}
		
		for (uint jj = 0; jj < active_pop_size; jj++) {
			j = active_idx[jj];
			if (positive_direction[j]) { // True
				if (cons_succ[j] >= DockConst_cons_limit) {
					rho[j] = LS_EXP_FACTOR * rho[j];
					cons_fail[j] = 0;
					cons_succ[j] = 0;
				}
				else if (cons_fail[j] >= DockConst_cons_limit) {
					rho[j] = LS_CONT_FACTOR * rho[j];
					cons_fail[j] = 0;
					cons_succ[j] = 0;
				}
				iteration_cnt[j]++;
			}
		}

		#if defined (PRINT_ALL)
		//printf("LS positive?: %u, iteration_cnt: %u, rho: %f, limit rho: %f\n",
		//positive_direction, iteration_cnt, rho, DockConst_rho_lower_bound);
		#endif
		// -----------------------------------------------

		float entity_possible_new_genotype[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];
		float genotype_bias[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];
		float deviate_plus_bias[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];
		float deviate_minus_bias[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];
		float rand_vec[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];  // = randf_vec...
		
		// Generating new random deviate
		// rho is the deviation of the uniform distribution
		for (uint i = 0; i < DockConst_num_of_genes; i++) {
			for (uint j = 0; j < pop_size; j++) {
				if (ls_is_active[j]) {
					float tmp_prng = rand_vec[j];

					// tmp1 is genotype_deviate
					float tmp1 = rho[j] * (2.0f * tmp_prng - 1.0f);

					if (i < 3) {
						tmp1 = tmp1 * DockConst_base_dmov_mul_sqrt3;
					}
					else {
						tmp1 = tmp1 * DockConst_base_dang_mul_sqrt3;
					}

					float deviate = 0.4f * tmp1;

					// tmp2 is the addition: genotype_deviate + genotype_bias
					// tmp3 is entity_possible_new_genotype
					float tmp_bias = (iteration_cnt[j] == 1)? 0.0f : genotype_bias[i][j];
					float bias = 0.6f * tmp_bias;

					deviate_plus_bias[i][j] = deviate + bias;
					deviate_minus_bias[i][j] = deviate - bias;

					float tmp2 = tmp1 + tmp_bias;
					float tmp3;
					if (positive_direction[j]) {
						tmp3 = genotype[i][j] + tmp2;
					} else {
						tmp3 = genotype[i][j] - tmp2;
					}

					if (i > 2) {
						if (i == 4) {
							tmp3 = map_angle_180(tmp3);
						}
						else {
							tmp3 = map_angle_360(tmp3);
						}
					}

					entity_possible_new_genotype[i][active_compr_idx[j]] = tmp3;
				}
			}
		}

		float energy_ia_ls[MAX_POPSIZE];
		float energy_ie_ls[MAX_POPSIZE];
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
			active_pop_size,
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
			active_pop_size,
			energy_ia_ls,
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
			active_pop_size,
			energy_ie_ls,
			local_coords_x,
			local_coords_y,
			local_coords_z
		);

		float candidate_energy[MAX_POPSIZE];
		for (uint jj = 0; jj < active_pop_size; jj++) {
			candidate_energy[jj] = energy_ia_ls[jj] + energy_ie_ls[jj];
		}

		// Updating LS energy-evaluation count
		LS_eval += active_pop_size;;

		for (uint jj = 0; jj < active_pop_size; jj++) {
			j = active_idx[jj];
			if (candidate_energy[jj] < current_energy[j]) {
				// Updating offspring_genotype & genotype_bias

				for (uint i = 0; i < DockConst_num_of_genes; i++) {
					genotype_bias[i][j] = positive_direction[j] ? deviate_plus_bias[i][j] : deviate_minus_bias[i][j];
					genotype[i][j] = entity_possible_new_genotype[i][jj];
				}

				current_energy[j] = candidate_energy[jj];
				cons_succ[j]++;
				cons_fail[j] = 0;
				positive_direction[j] = 1;
			}
			else {
				// Updating (halving) genotype_bias

				for (uint i = 0; i < DockConst_num_of_genes; i++) {
					genotype_bias[i][j] = (iteration_cnt[j] == 1)? 0.0f: (0.5f*genotype_bias[i][j]);
				}

				if (!positive_direction[j]) {
					cons_fail[j]++;
					cons_succ[j] = 0;
				}
				if (positive_direction[j]) {
					positive_direction[j] = 0;
				} else {
					positive_direction[j] = 1;
				}
			}
		}

		num_active_ls = pop_size;
		for (uint j = 0; j < pop_size; j++) {
			if ((iteration_cnt[j] > DockConst_max_num_of_iters) ||
			    (rho[j] <= DockConst_rho_lower_bound)) {
				ls_is_active[j] = 0;
				num_active_ls--;
			}
		}

	} // end of while (iteration_cnt) && (rho)

	// Writing resulting number of energy evals performed in LS
	*out_eval = LS_eval;

	// Writing resulting genotype and energy
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
		for (uint j = 0; j < pop_size; j++) {
			in_out_genotype[i] = genotype[i];
		}
	}
	for (uint j = 0; j < pop_size; j++) {
		in_out_energy[j] = current_energy[j];
	}
	
	#if defined (PRINT_ALL_IE) 
	printf("\n");
	printf("Finishing <local search>\n");
	printf("\n");
	#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
