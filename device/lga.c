#include "auxiliary.h"
#include "energy_ia.h"
#include "energy_ie.h"
#include "calc_pc.h"
#include "perform_ls.h"

/*
IC:  initial calculation of energy of populations
GG:  genetic generation
LS:  local search
*/

// --------------------------------------------------------------------------
// Lamarckian Genetic-Algorithm (GA): GA + LS (Local Search)
// --------------------------------------------------------------------------
void lga (
	const 	float*		PopulationCurrentInitial,
		float*  	PopulationCurrentFinal,
		float*  	EnergyCurrent,
		uint*		Evals_performed,
		uint*		Gens_performed,
		uint		DockConst_pop_size,
		uint            DockConst_num_of_energy_evals,
		uint            DockConst_num_of_generations,
		float           DockConst_tournament_rate,
		float           DockConst_mutation_rate,
		float           DockConst_abs_max_dmov,
		float           DockConst_abs_max_dang,
		float           Host_two_absmaxdmov,
		float           Host_two_absmaxdang,
		float           DockConst_crossover_rate,
		uchar           DockConst_num_of_genes,
	// pc
	const	int* 		PC_rotlist,
	const	float*		PC_ref_coords_x,// TODO: merge them into a single one?
	const	float*		PC_ref_coords_y,
	const	float*		PC_ref_coords_z,
	const	float*		PC_rotbonds_moving_vectors,
	const	float*		PC_rotbonds_unit_vectors,
	const	float*		PC_ref_orientation_quats,
			uint		DockConst_rotbondlist_length,
	// ia
	const 	float*		IA_IE_atom_charges,
	const	int*		IA_IE_atom_types,
	const	int*		IA_intraE_contributors,
	const	float*		IA_reqm,
	const	float*		IA_reqm_hbond,
	const	uint*		IA_atom1_types_reqm,
	const	uint*		IA_atom2_types_reqm,
	const	float*		IA_VWpars_AC,
	const	float*		IA_VWpars_BD,
	const	float*		IA_dspars_S,
	const	float*		IA_dspars_V,
			float		DockConst_smooth,
			uint		DockConst_num_of_intraE_contributors,
			float		DockConst_grid_spacing,
			uint		DockConst_num_of_atypes,
			float		DockConst_coeff_elec,
			float		DockConst_qasp,
			float		DockConst_coeff_desolv,
	// ie
	const	float*			Fgrids,
			uchar		DockConst_xsz,
			uint		DockConst_ysz,
			uint		DockConst_zsz,
			uchar		DockConst_num_of_atoms,
			float		DockConst_gridsize_x_minus1,
			float		DockConst_gridsize_y_minus1,
			float		DockConst_gridsize_z_minus1,
			uint		Host_mul_tmp2,
			uint		Host_mul_tmp3,
	// ls
			ushort		DockConst_max_num_of_iters,
			float		DockConst_rho_lower_bound,
			float		DockConst_base_dmov_mul_sqrt3,
			float		DockConst_base_dang_mul_sqrt3,
			uchar		DockConst_cons_limit,
	// Values changing every LGA run
			uint		Host_RunId,
			uint		Host_Offset_Pop,
			uint		Host_Offset_Ene
)
{
#if defined (PRINT_ALL_KRNL)
	printf("\n");
	printf("Starting <LGA run> ... \n");
	printf("\n");
	printf("%-40s %u\n", "DockConst_pop_size: ",        	DockConst_pop_size);
	printf("%-40s %u\n", "DockConst_num_of_energy_evals: ", DockConst_num_of_energy_evals);
	printf("%-40s %u\n", "DockConst_num_of_generations: ",  DockConst_num_of_generations);
	printf("%-40s %f\n", "DockConst_tournament_rate: ", 	DockConst_tournament_rate);
	printf("%-40s %f\n", "DockConst_mutation_rate: ", 		DockConst_mutation_rate);
	printf("%-40s +/-%fA\n", "DockConst_abs_max_dmov: ",	DockConst_abs_max_dmov);
	printf("%-40s +/-%f°\n", "DockConst_abs_max_dang: ",  	DockConst_abs_max_dang);
	printf("%-40s +/-%fA\n", "Host_two_absmaxdmov: ",		Host_two_absmaxdmov);
	printf("%-40s +/-%f°\n", "Host_two_absmaxdang: ",  		Host_two_absmaxdang);
	printf("%-40s %f\n", "DockConst_crossover_rate: ", 		DockConst_crossover_rate);
	printf("%-40s %u\n", "DockConst_num_of_genes: ",        DockConst_num_of_genes);
	printf("%-40s %u\n", "Host_RunId: ",        			Host_RunId);
	printf("%-40s %u\n", "Host_Offset_Pop: ",        		Host_Offset_Pop);
	printf("%-40s %u\n", "Host_Offset_Ene: ",        		Host_Offset_Ene);
	printf("\n");
	printf("%-40s %u\n", "DockConst_rotbondlist_length: ",  DockConst_rotbondlist_length);
	printf("\n");
	printf("%-40s %f\n", "DockConst_smooth: ", 				DockConst_smooth);
	printf("%-40s %u\n", "DockConst_num_of_intraE_contributors: ",	DockConst_num_of_intraE_contributors);
	printf("%-40s %f\n", "DockConst_grid_spacing: ", 		DockConst_grid_spacing);
	printf("%-40s %u\n", "DockConst_num_of_atypes: ",		DockConst_num_of_atypes);
	printf("%-40s %f\n", "DockConst_coeff_elec: ", 			DockConst_coeff_elec);
	printf("%-40s %f\n", "DockConst_qasp: ", 				DockConst_qasp);
	printf("%-40s %f\n", "DockConst_coeff_desolv: ", 		DockConst_coeff_desolv);
	printf("\n");
	printf("%-40s %u\n", "DockConst_xsz: ",					DockConst_xsz);
	printf("%-40s %u\n", "DockConst_ysz: ",					DockConst_ysz);
	printf("%-40s %u\n", "DockConst_zsz: ",					DockConst_zsz);
	printf("%-40s %u\n", "DockConst_num_of_atoms: ",		DockConst_num_of_atoms);
	printf("%-40s %f\n", "DockConst_gridsize_x_minus1: ", 	DockConst_gridsize_x_minus1);
	printf("%-40s %f\n", "DockConst_gridsize_y_minus1: ", 	DockConst_gridsize_y_minus1);
	printf("%-40s %f\n", "DockConst_gridsize_z_minus1: ", 	DockConst_gridsize_z_minus1);
	printf("%-40s %u\n", "Host_mul_tmp2: ",					Host_mul_tmp2);
	printf("%-40s %u\n", "Host_mul_tmp3: ",					Host_mul_tmp3);
	printf("\n");
	printf("%-40s %u\n", "DockConst_max_num_of_iters: ", 	DockConst_max_num_of_iters);
	printf("%-40s %f\n", "DockConst_rho_lower_bound: ", 	DockConst_rho_lower_bound);
	printf("%-40s %f\n", "DockConst_base_dmov_mul_sqrt3: ", DockConst_base_dmov_mul_sqrt3);
	printf("%-40s %f\n", "DockConst_base_dang_mul_sqrt3: ", DockConst_base_dang_mul_sqrt3);
	printf("%-40s %u\n", "DockConst_cons_limit: ",			DockConst_cons_limit);
#endif

	// --------------------------------------------------------------------------
	// ga
	float* GlobPopCurrFinal   = PopulationCurrentFinal;
	float* GlobEneCurr        = EnergyCurrent;
	uint* GlobEvals_performed = Evals_performed;
	uint* GlobGens_performed  = Gens_performed;
	// pc

	// ia

	// ie
	const	float*	IE_Fgrids = Fgrids;

	// --------------------------------------------------------------------------
	float LocalPopCurr[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];
	float LocalEneCurr[MAX_POPSIZE];

	float local_coords_x[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_y[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_z[MAX_NUM_OF_ATOMS][MAX_POPSIZE];

	for (uint i = 0; i < DockConst_num_of_atoms; i++) {
		for (int j = 0; j < MAX_POPSIZE; j++) {
			local_coords_x[i][j] = 0.0f;
			local_coords_y[i][j] = 0.0f;
			local_coords_z[i][j] = 0.0f;
		}
	}

	// --------------------------------------------------------------------------
	// Initial Calculation (IC) of scores
	// --------------------------------------------------------------------------
#if defined (PRINT_ALL_KRNL)
	printf("\n");
	printf("Starting <initial calculation> ... \n");
	printf("\n");
#endif

	// Read genotypes
	for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
		for (int pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
			float tmp_gene = PopulationCurrentInitial[Host_Offset_Pop + pop_cnt * ACTUAL_GENOTYPE_LENGTH + gene_cnt];
			LocalPopCurr[gene_cnt][pop_cnt] = tmp_gene;
		}
	}

	// Calculate energy
	float energy_ia_ic[MAX_POPSIZE];
	float energy_ie_ic[MAX_POPSIZE];
        
#if defined (PRINT_ALL_KRNL)
    printf("<before energy calc>\n");
#endif

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
		DockConst_pop_size,
		LocalPopCurr,
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
		//DockConst_num_of_atoms,
		DockConst_smooth,
		DockConst_num_of_intraE_contributors,
		DockConst_grid_spacing,
		DockConst_num_of_atypes,
		DockConst_coeff_elec,
		DockConst_qasp,
		DockConst_coeff_desolv,
		DockConst_pop_size,
		energy_ia_ic,
		local_coords_x,
		local_coords_y,
		local_coords_z
	);

	energy_ie(
		IE_Fgrids,
		IA_IE_atom_charges,
		IA_IE_atom_types,
		DockConst_xsz,
		DockConst_ysz,
		DockConst_zsz,
		DockConst_num_of_atoms,
		DockConst_gridsize_x_minus1,
		DockConst_gridsize_y_minus1,
		DockConst_gridsize_z_minus1,
		Host_mul_tmp2,
		Host_mul_tmp3,
		DockConst_pop_size,
		energy_ie_ic,
		local_coords_x,
		local_coords_y,
		local_coords_z
	);

	for (int pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
		LocalEneCurr[pop_cnt] = energy_ia_ic[pop_cnt] + energy_ie_ic[pop_cnt];
	}

#if defined (PRINT_ALL_KRNL)
	for (uint pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
		printf("Individual < after energy calc>: %3u, %20.6f\n", pop_cnt, LocalEneCurr[pop_cnt]);
	}
#endif

#if defined (PRINT_ALL_KRNL)
	printf("\n");
	printf("Finishing <initial calculation>\n");
	printf("\n");
#endif

	// --------------------------------------------------------------------------
	uint eval_cnt = DockConst_pop_size; // takes into account the IC evals
	uint generation_cnt = 0;

	while ((eval_cnt < DockConst_num_of_energy_evals) && (generation_cnt < DockConst_num_of_generations)) {

		float LocalPopNext[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];
		float LocalEneNext[MAX_POPSIZE];

		// --------------------------------------------------------------------------
		// Genetic Generation (GG)
		// --------------------------------------------------------------------------
#if defined (PRINT_ALL_KRNL)
		printf("\n");
		printf("Starting <finding best individual> ... \n");
		printf("\n");
#endif

		float loc_energies[MAX_POPSIZE];

		ushort best_entity = 0;
		loc_energies[0] = LocalEneCurr[0];

		for (int pop_cnt = 1; pop_cnt < DockConst_pop_size; pop_cnt++) {
			// copy energy to local memory
			loc_energies[pop_cnt] = LocalEneCurr[pop_cnt];

#if defined (PRINT_ALL_KRNL)
			if (pop_cnt == 1)  {
				printf("\n");
				printf("%3u %20.6f\n", 0, loc_energies[0]);
			}
			printf("%3u %20.6f\n", pop_cnt, loc_energies[pop_cnt]);
#endif

			if (loc_energies[pop_cnt] < loc_energies[best_entity]) {
				best_entity = pop_cnt;
			}
		}

#if defined (PRINT_ALL_KRNL)
		printf("best_entity: %3u, energy: %20.6f\n", best_entity, loc_energies[best_entity]);
#endif

#if defined (PRINT_ALL_KRNL)
		printf("\n");
		printf("Finishing <finding best individual>\n");
		printf("\n");
#endif

#if defined (PRINT_ALL_KRNL)
		printf("\n");
		printf("Starting <genetic generation> ... \n");
		printf("\n");
#endif

		float randv[10 + DockConst_num_of_genes][MAX_POPSIZE];
		randf_vec((float *)&randv[0][0], (10 + DockConst_num_of_genes) * MAX_POPSIZE);
                
		// ---------------------------------------------------
		// Elitism: copying the best entity to new population
		// ---------------------------------------------------
		for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
			LocalPopNext[gene_cnt][0] = LocalPopCurr[gene_cnt][best_entity];
		} 		
		LocalEneNext[0] = loc_energies[best_entity];

		for (int new_pop_cnt = 1; new_pop_cnt < DockConst_pop_size; new_pop_cnt++) {

			float local_entity_1 [ACTUAL_GENOTYPE_LENGTH];
			float local_entity_2 [ACTUAL_GENOTYPE_LENGTH]; 
		
			// ---------------------------------------------------
			// Binary-Tournament (BT) selection
			// ---------------------------------------------------

			// Get ushort binary_tournament selection prngs (parent index)
			uint bt_tmp_u0 = (uint) (DockConst_pop_size * randv[0][new_pop_cnt]);
			uint bt_tmp_u1 = (uint) (DockConst_pop_size * randv[1][new_pop_cnt]);
			uint bt_tmp_u2 = (uint) (DockConst_pop_size * randv[2][new_pop_cnt]);
			uint bt_tmp_u3 = (uint) (DockConst_pop_size * randv[3][new_pop_cnt]);

			// Get float binary_tournament selection prngs (tournament rate)
			float bt_tmp_f0 = randv[4][new_pop_cnt];
			float bt_tmp_f1 = randv[5][new_pop_cnt];
			float bt_tmp_f2 = randv[6][new_pop_cnt];
			float bt_tmp_f3 = randv[7][new_pop_cnt];

			uint parent1;
			uint parent2; 

			// First parent
			if (loc_energies[bt_tmp_u0] < loc_energies[bt_tmp_u1]) {
				if (bt_tmp_f0 < DockConst_tournament_rate) {
					parent1 = bt_tmp_u0;
				}
				else {
					parent1 = bt_tmp_u1;
				}
			}
			else {
				if (bt_tmp_f1 < DockConst_tournament_rate) {
					parent1 = bt_tmp_u1;
				}
				else {
					parent1 = bt_tmp_u0;
				}
			}

			// The better will be the second parent
			if (loc_energies[bt_tmp_u2] < loc_energies[bt_tmp_u3]) {
				if (bt_tmp_f2 < DockConst_tournament_rate) {
					parent2 = bt_tmp_u2;
				}
				else {
					parent2 = bt_tmp_u3;
				}
			}
			else {
				if (bt_tmp_f3 < DockConst_tournament_rate) {
					parent2 = bt_tmp_u3;
				}
				else {
					parent2 = bt_tmp_u2;
				}
			}

			// local_entity_1 and local_entity_2 are population-parent1, population-parent2
			for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
				local_entity_1[gene_cnt] = LocalPopCurr[gene_cnt][parent1];
				local_entity_2[gene_cnt] = LocalPopCurr[gene_cnt][parent2];
			}

			// --------------------------------------------------------------------------
			// Mating parents
			// --------------------------------------------------------------------------

			// get uchar genetic_generation prngs (gene index)
			// get float genetic_generation prngs (mutation rate)
			uint prng_GG_C_x = (uint) (DockConst_num_of_genes * randv[8][new_pop_cnt]);
			uint prng_GG_C_y = (uint) (DockConst_num_of_genes * randv[9][new_pop_cnt]);

			uint covr_point_low;
			uint covr_point_high;
			boolean twopoint_cross_yes = False;

			if (prng_GG_C_x == prng_GG_C_y) {
				covr_point_low = prng_GG_C_x;
			}
			else {
				twopoint_cross_yes = True;
				if (prng_GG_C_x < prng_GG_C_y) {
					covr_point_low = prng_GG_C_x;
					covr_point_high = prng_GG_C_y;
				}
				else {
					covr_point_low  = prng_GG_C_y;
					covr_point_high = prng_GG_C_x;
				}
			}
			
			// Reuse of bt prng float as crossover-rate
			boolean crossover_yes = False;
			if (DockConst_crossover_rate > bt_tmp_f0) {
				crossover_yes = True;
			}

			// Crossover and mutation
			for (int gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {

				// TODO: VERIFY IT
				float prngGG = randv[10 + gene_cnt][new_pop_cnt];
				float tmp_offspring;

				// Performing crossover
				if (   	(
					(crossover_yes == True) && (										// crossover
					( (twopoint_cross_yes == True)  && ((gene_cnt <= covr_point_low) || (gene_cnt > covr_point_high)) )  ||	// two-point crossover 			 		
					( (twopoint_cross_yes == False) && (gene_cnt <= covr_point_low))  					// one-point crossover
					)) || 
					(crossover_yes == False)	// no crossover
				   ) {
					tmp_offspring = local_entity_1[gene_cnt];
				}
				else {
					tmp_offspring = local_entity_2[gene_cnt];
				}

				// Performing mutation
				if (DockConst_mutation_rate > prngGG) {
					if(gene_cnt < 3) {
						tmp_offspring = tmp_offspring + Host_two_absmaxdmov * prngGG - DockConst_abs_max_dmov;
					}
					else {
						float tmp;
						tmp = tmp_offspring + Host_two_absmaxdang * prngGG - DockConst_abs_max_dang;
						if (gene_cnt == 4) {
							tmp_offspring = map_angle_180(tmp);
						}
						else {
							tmp_offspring = map_angle_360(tmp);
						}
					}
				}

				LocalPopNext[gene_cnt][new_pop_cnt] = tmp_offspring;
			}

        } // End loop over new_pop_cnt

        // Calculate energy
#if 1
        // make sure the buffers are 8 byte aligned
		double energy_gg[MAX_POPSIZE];
		float (*energy_ia_gg)[MAX_POPSIZE] = (void *)&energy_gg[0];
		float (*energy_ie_gg)[MAX_POPSIZE] = (void *)&energy_gg[MAX_POPSIZE/2];
#else
		float energy_ia_gg[MAX_POPSIZE];
		float energy_ie_gg[MAX_POPSIZE];
#endif

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
			DockConst_pop_size,
			LocalPopNext,     // Attention! in the original new_pop_cnt starts at 1, here it starts at 0!
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
			//DockConst_num_of_atoms,
			DockConst_smooth,
			DockConst_num_of_intraE_contributors,
			DockConst_grid_spacing,
			DockConst_num_of_atypes,
			DockConst_coeff_elec,
			DockConst_qasp,
			DockConst_coeff_desolv,
			DockConst_pop_size,     // Attention! in the original new_pop_cnt starts at 1, here it starts at 0!
			energy_ia_gg,
			local_coords_x,
			local_coords_y,
			local_coords_z
		);

		energy_ie(
			IE_Fgrids,
			IA_IE_atom_charges,
			IA_IE_atom_types,
			DockConst_xsz,
			DockConst_ysz,
			DockConst_zsz,
			DockConst_num_of_atoms,
			DockConst_gridsize_x_minus1,
			DockConst_gridsize_y_minus1,
			DockConst_gridsize_z_minus1,
			Host_mul_tmp2,
			Host_mul_tmp3,
			DockConst_pop_size,     // Attention! in the original new_pop_cnt starts at 1, here it starts at 0!
			energy_ie_gg,
			local_coords_x,
			local_coords_y,
			local_coords_z
		);

		for (uint new_pop_cnt = 0; new_pop_cnt < DockConst_pop_size; new_pop_cnt++) {
			LocalEneNext[new_pop_cnt] = (*energy_ia_gg)[new_pop_cnt] + (*energy_ie_gg)[new_pop_cnt];
        }

#if defined (PRINT_ALL_KRNL)
		for (uint new_pop_cnt = 0; new_pop_cnt < DockConst_pop_size; new_pop_cnt++) {
			printf("Individual < after energy calc>: %3u, %20.6f\n", new_pop_cnt, LocalEneNext[new_pop_cnt]);
		}
#endif
                
#if defined (PRINT_ALL_KRNL) 
		printf("\n");
		printf("Finishing <genetic generation>\n");
		printf("\n");
#endif

		// --------------------------------------------------------------------------
		// LS: Local Search
		// Performed on ALL individuals, not only for a small random subset
		// --------------------------------------------------------------------------

#if defined (PRINT_ALL_KRNL) 
		printf("\n");
		printf("Starting <local search> ... \n");
		printf("\n");
#endif

		uint ls_eval_cnt = 0;

#if defined (PRINT_ALL_KRNL)
		//printf("Individual <before ls>: %3u, %20.6f\n", entity_ls, LocalEneNext[entity_ls]);
#endif

		perform_ls(
			DockConst_max_num_of_iters,
			DockConst_rho_lower_bound,
			DockConst_base_dmov_mul_sqrt3,
			DockConst_num_of_genes,
			DockConst_base_dang_mul_sqrt3,
			DockConst_cons_limit,
			DockConst_pop_size,
			LocalPopNext,
			LocalEneNext,
			&ls_eval_cnt,

			PC_rotlist,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			DockConst_rotbondlist_length,
			Host_RunId,

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

			IE_Fgrids,
			DockConst_xsz,
			DockConst_ysz,
			DockConst_zsz,
			DockConst_num_of_atoms,
			DockConst_gridsize_x_minus1,
			DockConst_gridsize_y_minus1,
			DockConst_gridsize_z_minus1,
			Host_mul_tmp2,
			Host_mul_tmp3
		);

		// Accumulating number of evals
		//ls_eval_cnt += ls_eval_cnt_per_iter; // done inside perform_ls

#if defined (PRINT_ALL_KRNL)
		//printf("Individual < after ls>: %3u, %20.6f\n", entity_ls, LocalEneNext[entity_ls]);
#endif

#if defined (PRINT_ALL_KRNL)
		//printf("%u, ls_eval_cnt: %u\n", ls_ent_cnt, ls_eval_cnt);
#endif

#if defined (PRINT_ALL_KRNL)
		printf("\n");
		printf("Finishing <local search>\n");
		printf("\n");
#endif
		// ------------------------------------------------------------------

		// Update current pops & energies
		for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
			for (int pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
				LocalPopCurr[gene_cnt][pop_cnt] = LocalPopNext[gene_cnt][pop_cnt];
			}
		}

		for (int pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
			LocalEneCurr[pop_cnt] = LocalEneNext[pop_cnt];

#if defined (PRINT_ALL_KRNL)
			printf("Individual <updated>: %3u, %20.6f\n", pop_cnt, LocalEneCurr[pop_cnt]);
#endif
		}

		// Update energy evaluations count: count LS and GG evals
		eval_cnt += ls_eval_cnt + DockConst_pop_size; 

		// Update generation count
		generation_cnt++;

#if defined (PRINT_ALL_KRNL)
		printf("eval_cnt: %u, generation_cnt: %u\n", eval_cnt, generation_cnt);
#endif

	} // End while eval_cnt & generation_cnt

	// --------------------------------------------------------------------------
	// Write final pop & energies back to SX-Aurora VE external memory
	// --------------------------------------------------------------------------
	for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
		//printf("\n");
		//printf("pop_cnt: %u\n", pop_cnt);
		for (int pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
			GlobPopCurrFinal[Host_Offset_Pop + pop_cnt*ACTUAL_GENOTYPE_LENGTH + gene_cnt] = LocalPopCurr[gene_cnt][pop_cnt];
			//printf("Final geno: %3u, %20.6f\n", gene_cnt, LocalPopCurr[pop_cnt][gene_cnt]);
		}
    }

	for (int pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
		GlobEneCurr[Host_Offset_Ene + pop_cnt] = LocalEneCurr[pop_cnt];
		//printf("Final energy: %3u, %20.6f\n", pop_cnt, LocalEneCurr[pop_cnt]);
	}

	// Write final evals & generation counts to SX-Aurora VE external memory
	GlobEvals_performed[Host_RunId] = eval_cnt;
	GlobGens_performed [Host_RunId] = generation_cnt;

#if defined (PRINT_ALL_KRNL)
	printf("Host_RunId: %u, eval_cnt: %u, generation_cnt: %u\n", Host_RunId, eval_cnt, generation_cnt);
#endif

#if defined (PRINT_ALL_KRNL)
	printf("\n");
	printf("Finishing <LGA run>\n");
	printf("\n");
#endif
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
