#include "auxiliary.c"
#include "energy_ia.c"
#include "energy_ie.c"
#include "calc_pc.c"
#include "perform_ls.c"

/*
IC:  initial calculation of energy of populations
GG:  genetic generation
LS:  local search
*/

// --------------------------------------------------------------------------
// Lamarckian Genetic-Algorithm (GA): GA + LS (Local Search)
// --------------------------------------------------------------------------
void lga (
	const 	uint64_t	VEVMA_PopulationCurrentInitial,
			uint64_t  	VEVMA_PopulationCurrentFinal,
			uint64_t  	VEVMA_EnergyCurrent,
			uint64_t	VEVMA_Evals_performed,
            uint64_t	VEVMA_Gens_performed,
			uint64_t    VEVMA_dockpars_prng_states,
			uint		DockConst_pop_size,
		    uint        DockConst_num_of_energy_evals,
			uint        DockConst_num_of_generations,
		    float       DockConst_tournament_rate,
			float       DockConst_mutation_rate,
		    float       DockConst_abs_max_dmov,
			float       DockConst_abs_max_dang,
		    float       Host_two_absmaxdmov,
			float       Host_two_absmaxdang,
			float       DockConst_crossover_rate,
			uint        DockConst_num_of_lsentities,
			uchar       DockConst_num_of_genes,
	// pc
	const	uint64_t	VEVMA_pc_rotlist,
	const	uint64_t	VEVMA_pc_ref_coords_x,// TODO: merge them into a single one?
	const	uint64_t	VEVMA_pc_ref_coords_y,
	const	uint64_t	VEVMA_pc_ref_coords_z,
	const	uint64_t	VEVMA_pc_rotbonds_moving_vectors,
	const	uint64_t	VEVMA_pc_rotbonds_unit_vectors,
	const	uint64_t	VEVMA_pc_ref_orientation_quats,
			uint		DockConst_rotbondlist_length,
	// ia
	const 	uint64_t	VEVMA_ia_ie_atom_charges,
	const	uint64_t	VEVMA_ia_ie_atom_types,
	const	uint64_t	VEVMA_ia_intraE_contributors,
	const	uint64_t	VEVMA_ia_reqm,
	const	uint64_t	VEVMA_ia_reqm_hbond,
	const	uint64_t	VEVMA_ia_atom1_types_reqm,
	const	uint64_t	VEVMA_ia_atom2_types_reqm,
	const	uint64_t	VEVMA_ia_VWpars_AC,
	const	uint64_t	VEVMA_ia_VWpars_BD,
	const	uint64_t	VEVMA_ia_dspars_S,
	const	uint64_t	VEVMA_ia_dspars_V,
			float		DockConst_smooth,
			uint		DockConst_num_of_intraE_contributors,
			float		DockConst_grid_spacing,
			uint		DockConst_num_of_atypes,
			float		DockConst_coeff_elec,
			float		DockConst_qasp,
			float		DockConst_coeff_desolv,
	// ie
	const	uint64_t	VEVMA_Fgrids,
			uchar		DockConst_g1,
			uint		DockConst_g2,
			uint		DockConst_g3,
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
			uint 	    Host_Offset_Pop,
			uint	    Host_Offset_Ene
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
	printf("%-40s %u\n", "DockConst_num_of_lsentities: ",   DockConst_num_of_lsentities);
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
	printf("%-40s %u\n", "DockConst_g1: ",					DockConst_g1);
	printf("%-40s %u\n", "DockConst_g2: ",					DockConst_g2);
	printf("%-40s %u\n", "DockConst_g3: ",					DockConst_g3);
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
	const float* GlobPopCurrInitial = (float*)(VEVMA_PopulationCurrentInitial);
	      float* GlobPopCurrFinal   = (float*)(VEVMA_PopulationCurrentFinal);
	      float* GlobEneCurr        = (float*)(VEVMA_EnergyCurrent);
		  uint* GlobEvals_performed = (uint*)(VEVMA_Evals_performed);
		  uint* GlobGens_performed  = (uint*)(VEVMA_Gens_performed);
		  uint* dockpars_prng_states = (uint*) VEVMA_dockpars_prng_states;
	// pc
	const	int*	PC_rotlist = (int*)(VEVMA_pc_rotlist);
	const	float*	PC_ref_coords_x = (float*)(VEVMA_pc_ref_coords_x);
	const	float*	PC_ref_coords_y = (float*)(VEVMA_pc_ref_coords_y);
	const	float*	PC_ref_coords_z = (float*)(VEVMA_pc_ref_coords_z);
	const	float*	PC_rotbonds_moving_vectors = (float*)(VEVMA_pc_rotbonds_moving_vectors);
	const	float*	PC_rotbonds_unit_vectors = (float*)(VEVMA_pc_rotbonds_unit_vectors);
	const	float*	PC_ref_orientation_quats = (float*)(VEVMA_pc_ref_orientation_quats);
	// ia
	const	float*	IA_IE_atom_charges = (float*)(VEVMA_ia_ie_atom_charges);
	const	int*	IA_IE_atom_types = (int*)(VEVMA_ia_ie_atom_types);
	const	int*	IA_intraE_contributors = (int*)(VEVMA_ia_intraE_contributors);
	const	float*	IA_reqm = (float*)(VEVMA_ia_reqm);
	const	float*	IA_reqm_hbond = (float*)(VEVMA_ia_reqm_hbond);
	const	uint*	IA_atom1_types_reqm = (uint*)(VEVMA_ia_atom1_types_reqm);
	const	uint*	IA_atom2_types_reqm = (uint*)(VEVMA_ia_atom2_types_reqm);
	const	float*	IA_VWpars_AC = (float*)(VEVMA_ia_VWpars_AC);
	const	float*	IA_VWpars_BD = (float*)(VEVMA_ia_VWpars_BD);
	const	float*	IA_dspars_S = (float*)(VEVMA_ia_dspars_S);
	const	float*	IA_dspars_V = (float*)(VEVMA_ia_dspars_V);
	// ie
	const	float*	IE_Fgrids = (float*)(VEVMA_Fgrids);

	// --------------------------------------------------------------------------
	float LocalPopCurr[MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];
	float LocalEneCurr[MAX_POPSIZE];

	float local_coords_x[MAX_NUM_OF_ATOMS];
	float local_coords_y[MAX_NUM_OF_ATOMS];
	float local_coords_z[MAX_NUM_OF_ATOMS];

	uint prng_offset = Host_Offset_Ene; // run_cnt * dockpars.pop_siz

	// --------------------------------------------------------------------------
	// Initial Calculation (IC) of scores
	// --------------------------------------------------------------------------
	#if defined (PRINT_ALL_KRNL) 
	printf("\n");
	printf("Starting <initial calculation> ... \n");
	printf("\n");
	#endif

	for (uint pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {

		// Read genotype
		for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
			float tmp_gene = GlobPopCurrInitial[Host_Offset_Pop + pop_cnt * ACTUAL_GENOTYPE_LENGTH + gene_cnt];
			LocalPopCurr[pop_cnt][gene_cnt] = tmp_gene;
		}

		#if defined (PRINT_ALL_KRNL)
		printf("Individual <before energy calc>: %3u\n", pop_cnt);
		#endif

		// Calculate energy
		float energy_ia_ic;
		float energy_ie_ic;
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
			LocalPopCurr[pop_cnt],
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
			&energy_ia_ic,
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
			&energy_ie_ic,
			local_coords_x,
			local_coords_y,
			local_coords_z
		);
		LocalEneCurr[pop_cnt] = energy_ia_ic + energy_ie_ic;

		#if defined (PRINT_ALL_KRNL)
		printf("Individual < after energy calc>: %3u, %20.6f\n", pop_cnt, LocalEneCurr[pop_cnt]);
		#endif
	}

	#if defined (PRINT_ALL_KRNL)
	printf("\n");
	printf("Finishing <initial calculation>\n");
	printf("\n");
	#endif

	// --------------------------------------------------------------------------
	uint eval_cnt = DockConst_pop_size; // takes into account the IC evals
	uint generation_cnt = 0;

	while ((eval_cnt < DockConst_num_of_energy_evals) && (generation_cnt < DockConst_num_of_generations)) {

		float LocalPopNext[MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];
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

		for (uint pop_cnt = 1; pop_cnt < DockConst_pop_size; pop_cnt++) {
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

		/*
		#pragma ivdep array (LocalPopNext)
		#pragma ivdep array (LocalEneNext)
		*/
		for (uint new_pop_cnt = 1; new_pop_cnt < DockConst_pop_size; new_pop_cnt++) {

			// ---------------------------------------------------
			// Elitism: copying the best entity to new population
			// ---------------------------------------------------
			if (new_pop_cnt == 1) {
				for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
					LocalPopNext[0][gene_cnt] = LocalPopCurr[best_entity][gene_cnt];
				} 		
				LocalEneNext[0] = loc_energies[best_entity];
			}

			#if defined (PRINT_ALL_KRNL)
			#endif

			float local_entity_1 [ACTUAL_GENOTYPE_LENGTH];
			float local_entity_2 [ACTUAL_GENOTYPE_LENGTH]; 
		
			// ---------------------------------------------------
			// Binary-Tournament (BT) selection
			// ---------------------------------------------------

			// TODO: FIX INDEXES
			// TODO: VECTORIZE IT
			// Get ushort binary_tournament selection prngs (parent index)
			ushort bt_tmp_u0 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[prng_offset + new_pop_cnt]));
			ushort bt_tmp_u1 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[prng_offset + new_pop_cnt]));
			ushort bt_tmp_u2 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[prng_offset + new_pop_cnt]));
			ushort bt_tmp_u3 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[prng_offset + new_pop_cnt]));

			// TODO: FIX INDEXES
			// TODO: VECTORIZE IT
			// Get float binary_tournament selection prngs (tournament rate)
			float bt_tmp_f0 = randf(&dockpars_prng_states[prng_offset + new_pop_cnt]);
			float bt_tmp_f1 = randf(&dockpars_prng_states[prng_offset + new_pop_cnt]);
			float bt_tmp_f2 = randf(&dockpars_prng_states[prng_offset + new_pop_cnt]);
			float bt_tmp_f3 = randf(&dockpars_prng_states[prng_offset + new_pop_cnt]);

			ushort parent1;
			ushort parent2; 

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
				local_entity_1[gene_cnt] = LocalPopCurr[parent1][gene_cnt];
				local_entity_2[gene_cnt] = LocalPopCurr[parent2][gene_cnt];
			}

			// --------------------------------------------------------------------------
			// Mating parents
			// --------------------------------------------------------------------------

			// TODO: VECTORIZE IT
			// get uchar genetic_generation prngs (gene index)
			// get float genetic_generation prngs (mutation rate)
			uchar prng_GG_C_x = (uchar) (DockConst_num_of_genes * randf(&dockpars_prng_states[prng_offset + new_pop_cnt]));
			uchar prng_GG_C_y = (uchar) (DockConst_num_of_genes * randf(&dockpars_prng_states[prng_offset + new_pop_cnt]));

			uchar covr_point_low;
			uchar covr_point_high;
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
			for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {

				// TODO: VERIFY IT
				float prngGG = randf(&dockpars_prng_states[prng_offset + new_pop_cnt]);
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

				LocalPopNext[new_pop_cnt][gene_cnt] = tmp_offspring;
			}

			#if defined (PRINT_ALL_KRNL)	
			#endif	

			// Calculate energy
			float energy_ia_gg;
			float energy_ie_gg;
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
				LocalPopNext[new_pop_cnt],
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
				&energy_ia_gg,
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
				&energy_ie_gg,
				local_coords_x,
				local_coords_y,
				local_coords_z
			);
			LocalEneNext[new_pop_cnt] = energy_ia_gg + energy_ie_gg;

			#if defined (PRINT_ALL_KRNL)
			printf("Individual < after energy calc>: %3u, %20.6f\n", new_pop_cnt, LocalEneNext[new_pop_cnt]);
			#endif
		} 

		#if defined (PRINT_ALL_KRNL) 
		printf("\n");
		printf("Finishing <genetic generation>\n");
		printf("\n");
		#endif

		// --------------------------------------------------------------------------
		// LS: Local Search
		// Subject num_of_entity_for_ls pieces of offsprings to LS 
		// --------------------------------------------------------------------------

		#if defined (PRINT_ALL_KRNL) 
		printf("\n");
		printf("Starting <local search> ... \n");
		printf("\n");
		#endif

		uint ls_eval_cnt = 0;

		/*
		#pragma ivdep
		*/
		for (uint ls_ent_cnt = 0; ls_ent_cnt < DockConst_num_of_lsentities; ls_ent_cnt++) {

			uint ls_eval_cnt_per_iter;

			// TODO: FIX INDEX FOR PRNG
			// Choose random & different entities on every iteration
			ushort entity_ls = (ushort)(DockConst_pop_size * randf(&dockpars_prng_states[prng_offset]));

			#if defined (PRINT_ALL_KRNL)
			printf("Individual <before ls>: %3u, %20.6f\n", entity_ls, LocalEneNext[entity_ls]);
			#endif

			perform_ls(
				DockConst_max_num_of_iters,
				DockConst_rho_lower_bound,
				DockConst_base_dmov_mul_sqrt3,
				DockConst_num_of_genes,
				DockConst_base_dang_mul_sqrt3,
				DockConst_cons_limit,

				LocalPopNext[entity_ls],
				&LocalEneNext[entity_ls],
				&ls_eval_cnt_per_iter,
				&dockpars_prng_states[prng_offset + entity_ls],

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
				DockConst_g1,
				DockConst_g2,
				DockConst_g3,
				DockConst_num_of_atoms,
				DockConst_gridsize_x_minus1,
				DockConst_gridsize_y_minus1,
				DockConst_gridsize_z_minus1,
				Host_mul_tmp2,
				Host_mul_tmp3
			);

			// Accumulating number of evals
			ls_eval_cnt += ls_eval_cnt_per_iter;

			#if defined (PRINT_ALL_KRNL)
			printf("Individual < after ls>: %3u, %20.6f\n", entity_ls, LocalEneNext[entity_ls]);
			#endif

			#if defined (PRINT_ALL_KRNL)
			printf("%u, ls_eval_cnt: %u\n", ls_ent_cnt, ls_eval_cnt);
			#endif
		} // End of for-loop ls_ent_cnt

		#if defined (PRINT_ALL_KRNL) 
		printf("\n");
		printf("Finishing <local search>\n");
		printf("\n");
		#endif
		// ------------------------------------------------------------------

		// Update current pops & energies
		for (uint pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
			for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
				LocalPopCurr[pop_cnt][gene_cnt] = LocalPopNext[pop_cnt][gene_cnt];
			}

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
	// Write final pop & energies back to FPGA-board DDRs
	// --------------------------------------------------------------------------
	for (uint pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
		//printf("\n");
		//printf("pop_cnt: %u\n", pop_cnt);
		for (uint gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
			GlobPopCurrFinal[Host_Offset_Pop + pop_cnt*ACTUAL_GENOTYPE_LENGTH + gene_cnt] = LocalPopCurr[pop_cnt][gene_cnt];

			//printf("Final geno: %3u, %20.6f\n", gene_cnt, LocalPopCurr[pop_cnt][gene_cnt]);
		}

		GlobEneCurr[Host_Offset_Ene + pop_cnt] = LocalEneCurr[pop_cnt];
		//printf("Final energy: %3u, %20.6f\n", pop_cnt, LocalEneCurr[pop_cnt]);
	}

	// Write final evals & generation counts to FPGA-board DDRs
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
