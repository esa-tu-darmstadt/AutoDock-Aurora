#include <stdio.h>
#include <stdint.h>

#include "auxiliary.c"
#include "energy_ia.c"
#include "energy_ie.c"
#include "calc_pc.c"

#include "math.h"

//IC:  initial calculation of energy of populations
//GG:  genetic generation 
//LS:  local search
//OFF: turn off 

/*
 * -----------------------------------------------
 * Lamarckian Genetic-Algorithm (GA): GA + LS (Local Search)
 * -----------------------------------------------
 * */
uint64_t libkernel_ga (
	const 	uint64_t	VEVMA_PopulationCurrentInitial,
			uint64_t  	VEVMA_PopulationCurrentFinal,
			uint64_t  	VEVMA_EnergyCurrent,
			uint64_t	VEVMA_Evals_performed,
            uint64_t	VEVMA_Gens_performed,
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
	        ushort      Host_RunId,
			uint 	    Host_Offset_Pop,
			uint	    Host_Offset_Ene,
			uint64_t    VEVMA_dockpars_prng_states,
	/*
	 * pc
	 */
	const	uint64_t	VEVMA_Static_rotlist,
	const	uint64_t	VEVMA_Static_ref_coords_x,
	const	uint64_t	VEVMA_Static_ref_coords_y,
	const	uint64_t	VEVMA_Static_ref_coords_z,
	const	uint64_t	VEVMA_Static_rotbonds_moving_vectors,
	const	uint64_t	VEVMA_Static_rotbonds_unit_vectors,
			uint		DockConst_rotbondlist_length,
	//		uchar				DockConst_num_of_genes,
	const	uint64_t	VEVMA_Static_ref_orientation_quats,
	//		ushort				Host_RunId,
)
{
	#ifdef PRINT_ALL_KRNL
  	printf("Starting libkernel_ga ... \n");
  	#endif



	#if defined (DEBUG_KRNL_GA)
	printf("\n");
	printf("%-40s %u\n", "DockConst_pop_size: ",        		DockConst_pop_size);
	printf("%-40s %u\n", "DockConst_num_of_energy_evals: ",  	DockConst_num_of_energy_evals);
	printf("%-40s %u\n", "DockConst_num_of_generations: ",  	DockConst_num_of_generations);
	printf("%-40s %f\n", "DockConst_tournament_rate: ", 		DockConst_tournament_rate);
	printf("%-40s %f\n", "DockConst_mutation_rate: ", 		DockConst_mutation_rate);
	printf("%-40s +/-%fA\n", "DockConst_abs_max_dmov: ",		DockConst_abs_max_dmov);
	printf("%-40s +/-%f°\n", "DockConst_abs_max_dang: ",  		DockConst_abs_max_dang);
	printf("%-40s +/-%fA\n", "Host_two_absmaxdmov: ",		Host_two_absmaxdmov);
	printf("%-40s +/-%f°\n", "Host_two_absmaxdang: ",  		Host_two_absmaxdang);
	printf("%-40s %f\n", "DockConst_crossover_rate: ", 		DockConst_crossover_rate);
	printf("%-40s %u\n", "DockConst_num_of_lsentities: ",   	DockConst_num_of_lsentities);
	printf("%-40s %u\n", "DockConst_num_of_genes: ",        	DockConst_num_of_genes);
	#endif

	// ------------------------------------------------------------------

	/*
	 * ga
	 */
	const float* GlobPopCurrInitial = (float*)(VEVMA_PopulationCurrentInitial + Host_Offset_Pop);
	      float* GlobPopCurrFinal   = (float*)(VEVMA_PopulationCurrentFinal + Host_Offset_Pop);
	      float* GlobEneCurr        = (float*)(VEVMA_EnergyCurrent + Host_Offset_Ene);
		  uint* GlobEvals_performed = (uint*)(VEVMA_Evals_performed);
		  uint* GlobGens_performed  = (uint*)(VEVMA_Gens_performed);

	uint* dockpars_prng_states = (uint*) VEVMA_dockpars_prng_states;

	/*
	 * pc
	 */
	const	int*	KerConstStatic_rotlist = (int*)(VEVMA_Static_rotlist);
	const	float*	KerConstStatic_ref_coords_x = (float*)(VEVMA_Static_ref_coords_x);
	const	float*	KerConstStatic_ref_coords_y = (float*)(VEVMA_Static_ref_coords_y);
	const	float*	KerConstStatic_ref_coords_z = (float*)(VEVMA_Static_ref_coords_z);
	const	float*	KerConstStatic_rotbonds_moving_vectors = (float*)(VEVMA_Static_rotbonds_moving_vectors);
	const	float*	KerConstStatic_rotbonds_unit_vectors = (float*)(VEVMA_Static_rotbonds_unit_vectors);
	const	float*	KerConstStatic_ref_orientation_quats = (float*)(VEVMA_Static_ref_orientation_quats);

	// ------------------------------------------------------------------

	float LocalPopCurr[MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];
	float LocalEneCurr[MAX_POPSIZE];

	// ------------------------------------------------------------------
	// Initial Calculation (IC) of scores
	// ------------------------------------------------------------------
	for (ushort pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
		// Calculate energy
		for (uchar gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
			float tmp_gene = GlobPopCurrInitial[pop_cnt*ACTUAL_GENOTYPE_LENGTH + gene_cnt];
			LocalPopCurr[pop_cnt][gene_cnt] = tmp_gene;
		}

		float energy_ia_ic;
		float energy_ie_ic;
		calc_pc(
			KerConstStatic_rotlist,
			KerConstStatic_ref_coords_x,
			KerConstStatic_ref_coords_y,
			KerConstStatic_ref_coords_z,
			KerConstStatic_rotbonds_moving_vectors,
			KerConstStatic_rotbonds_unit_vectors,
			DockConst_rotbondlist_length,
			DockConst_num_of_genes,
			KerConstStatic_ref_orientation_quats,
			Host_RunId,
			LocalPopCurr[pop_cnt],

		);
		energy_ia();
		energy_ie();
		LocalEneCurr[pop_cnt] = energy_ia_ic + energy_ie_ic;
	}
	// ------------------------------------------------------------------

	uint eval_cnt = DockConst_pop_size; // takes into account the IC evals

	uint generation_cnt = 0;

	while ((eval_cnt < DockConst_num_of_energy_evals) && (generation_cnt < DockConst_num_of_generations)) {

		float LocalPopNext[MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];
		float LocalEneNext[MAX_POPSIZE];

		// ------------------------------------------------------------------
		// Genetic Generation (GG)
		// ------------------------------------------------------------------
		float loc_energies[MAX_POPSIZE];
		ushort best_entity = 0;

		for (ushort pop_cnt = 1; pop_cnt < DockConst_pop_size; pop_cnt++) {
			// copy energy to local memory
			loc_energies[pop_cnt] = LocalEneCurr[pop_cnt];

			#if defined (DEBUG_KRNL_GA)
			if (pop_cnt==0) {printf("\n");}
			printf("%3u %20.6f\n", pop_cnt, loc_energies[pop_cnt]);
			#endif

			if (loc_energies[pop_cnt] < loc_energies[best_entity]) {
				best_entity = pop_cnt;
			}
		}

		#if defined (DEBUG_KRNL_GA)
		printf("best_entity: %3u, energy: %20.6f\n", best_entity, loc_energies[best_entity]);
		#endif

		/*
		#pragma ivdep array (LocalPopNext)
		#pragma ivdep array (LocalEneNext)
		*/
		for (ushort new_pop_cnt = 1; new_pop_cnt < DockConst_pop_size; new_pop_cnt++) {

			// ---------------------------------------------------
			// Elitism: copying the best entity to new population
			// ---------------------------------------------------
			if (new_pop_cnt == 1) {
				for (uchar gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
					LocalPopNext[0][gene_cnt] = LocalPopCurr[best_entity][gene_cnt];
				} 		
				LocalEneNext[0] = loc_energies[best_entity];
			}

			#if defined (DEBUG_KRNL_GA)
			printf("Krnl_GA: %u\n", new_pop_cnt);
			#endif

			float local_entity_1 [ACTUAL_GENOTYPE_LENGTH];
			float local_entity_2 [ACTUAL_GENOTYPE_LENGTH]; 
		
			// ---------------------------------------------------
			// Binary-Tournament (BT) selection
			// ---------------------------------------------------

			// TODO: FIX INDEXES
			// TODO: VECTORIZE IT
			// Get ushort binary_tournament selection prngs (parent index)
			ushort bt_tmp_u0 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[new_pop_cnt]));
			ushort bt_tmp_u1 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[new_pop_cnt]));
			ushort bt_tmp_u2 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[new_pop_cnt + 1]));
			ushort bt_tmp_u3 = (ushort) (DockConst_pop_size * randf(&dockpars_prng_states[new_pop_cnt + 1]));

			// TODO: FIX INDEXES
			// TODO: VECTORIZE IT
			// Get float binary_tournament selection prngs (tournament rate)
			float bt_tmp_f0 = randf(&dockpars_prng_states[new_pop_cnt + 2]);
			float bt_tmp_f1 = randf(&dockpars_prng_states[new_pop_cnt + 2]);
			float bt_tmp_f2 = randf(&dockpars_prng_states[new_pop_cnt + 3]);
			float bt_tmp_f3 = randf(&dockpars_prng_states[new_pop_cnt + 3]);

			ushort parent1;
			ushort parent2; 

			// First parent
			if (loc_energies[bt_tmp_u0] < loc_energies[bt_tmp_u1]) {
				if (bt_tmp_f0 < DockConst_tournament_rate) {parent1 = bt_tmp_u0;}
				else				           {parent1 = bt_tmp_u1;}
			}
			else {
				if (bt_tmp_f1 < DockConst_tournament_rate) {parent1 = bt_tmp_u1;}
				else				           {parent1 = bt_tmp_u0;}
			}

			// The better will be the second parent
			if (loc_energies[bt_tmp_u2] < loc_energies[bt_tmp_u3]) {
				if (bt_tmp_f2 < DockConst_tournament_rate) {parent2 = bt_tmp_u2;}
				else		          	           {parent2 = bt_tmp_u3;}
			}
			else {
				if (bt_tmp_f3 < DockConst_tournament_rate) {parent2 = bt_tmp_u3;}
				else			                   {parent2 = bt_tmp_u2;}
			}

			// local_entity_1 and local_entity_2 are population-parent1, population-parent2
			for (uchar gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
				local_entity_1[gene_cnt] = LocalPopCurr[parent1][gene_cnt];
				local_entity_2[gene_cnt] = LocalPopCurr[parent2][gene_cnt];
			}

			// ---------------------------------------------------
			// Mating parents
			// ---------------------------------------------------	

			// TODO: VECTORIZE IT
			// get uchar genetic_generation prngs (gene index)
			// get float genetic_generation prngs (mutation rate)
			uchar prng_GG_C_x = (uchar) (DockConst_num_of_genes * randf(&dockpars_prng_states[new_pop_cnt + 4]));
			uchar prng_GG_C_y = (uchar) (DockConst_num_of_genes * randf(&dockpars_prng_states[new_pop_cnt + 4]));

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
			boolean crossover_yes = (DockConst_crossover_rate > bt_tmp_f0);

			// Crossover and mutation
			for (uchar gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {

				// TODO: VERIFY IT
				float prngGG = randf(&dockpars_prng_states[new_pop_cnt + 5]);
				float tmp_offspring;

				// Performing crossover
				if (   	(
					crossover_yes && (										// crossover
					( (twopoint_cross_yes == True)  && ((gene_cnt <= covr_point_low) || (gene_cnt > covr_point_high)) )  ||	// two-point crossover 			 		
					( (twopoint_cross_yes == False) && (gene_cnt <= covr_point_low))  					// one-point crossover
					)) || 
					(!crossover_yes)	// no crossover
				   ) {
					tmp_offspring = local_entity_1[gene_cnt];
				}
				else {
					tmp_offspring = local_entity_2[gene_cnt];
				}

				// Performing mutation
				if (DockConst_mutation_rate > prngGG) {
					if(gene_cnt < 3) {
						tmp_offspring = tmp_offspring + Host_two_absmaxdmov*prngGG-DockConst_abs_max_dmov;
					}
					else {
						float tmp;
						tmp = tmp_offspring + Host_two_absmaxdang*prngGG-DockConst_abs_max_dang;
						if (gene_cnt == 4) {
							tmp_offspring = map_angle_180(tmp);
						}
						else {
							tmp_offspring = map_angle_360(tmp);
						}
					}
				}

				// Calculate energy
				LocalPopNext [new_pop_cnt][gene_cnt] = tmp_offspring;
			}

			#if defined (DEBUG_KRNL_GG)
			printf("GG - tx pop: %u", new_pop_cnt); 		
			#endif	

			// Read energy
			float energy_ia_gg;
			float energy_ie_gg;
			// TODO: CALC ENERGY
			LocalEneNext[new_pop_cnt] = energy_ia_gg + energy_ie_gg;

			#if defined (DEBUG_KRNL_GG)
			printf(", GG - rx pop: %u\n", new_pop_cnt); 		
			#endif
		} 
		// ------------------------------------------------------------------
		// LS: Local Search
		// Subject num_of_entity_for_ls pieces of offsprings to LS 
		// ------------------------------------------------------------------

		uint ls_eval_cnt = 0;

		/*
		#pragma ivdep
		*/
		for (ushort ls_ent_cnt = 0; ls_ent_cnt < DockConst_num_of_lsentities; ls_ent_cnt++) {

			// TODO: FIX INDEX FOR PRNG
			// Choose random & different entities on every iteration
			ushort entity_ls = (ushort)(DockConst_num_of_lsentities * randf(&dockpars_prng_states[ls_ent_cnt]));

/*			
			write_pipe_block(pipe00ga2ls00ls100energy, &LocalEneNext[entity_ls1]);
			// LOOP_GA_LS_INNER_WRITE_GENOTYPE
			for (uchar gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
				write_pipe_block(pipe00ga2ls00ls100genotype, &LocalPopNext[entity_ls1][gene_cnt]);
			}
*/

			// TODO: FIX THESE DECLARATIONS
/*
			float2 evalenergy_tmp1;
*/

		
			#if defined (DEBUG_KRNL_LS)
			printf("LS - got all eval & energies back\n");
			#endif

			// TODO: RETURNING ENERGYCALC COUNT
/*
			float eetmp1 = evalenergy_tmp1.x;
			uint eval_tmp1 = *(uint*)&eetmp1;
*/
			// TODO: RETURNING ENERGIES CALCULATED
/*			
			LocalEneNext[entity_ls1] = evalenergy_tmp1.y;
*/
			// TODO: READ RETURNING GENOTYPES
/*			
			for (uchar gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
				read_pipe_block(pipe00ls2ga00ls100genotype, &LocalPopNext[entity_ls1][gene_cnt]);
			}
*/
			// TODO: SUM UP EVALS
/*
			ls_eval_cnt += eval_tmp1 + eval_tmp2 + eval_tmp3 + eval_tmp4 + eval_tmp5 + eval_tmp6 + eval_tmp7 + eval_tmp8 + eval_tmp9;
*/
			#if defined (DEBUG_KRNL_LS)
			printf("%u, ls_eval_cnt: %u\n", ls_ent_cnt, ls_eval_cnt);
			printf("LS - got all genotypes back\n");
			#endif
		} // End of for-loop ls_ent_cnt
		// ------------------------------------------------------------------

		// Update current pops & energies
		for (ushort pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
			for (uchar gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
				LocalPopCurr[pop_cnt][gene_cnt] = LocalPopNext[pop_cnt][gene_cnt];
			}

			LocalEneCurr[pop_cnt] = LocalEneNext[pop_cnt];
		}

		// Update energy evaluations count: count LS and GG evals
		eval_cnt += ls_eval_cnt + DockConst_pop_size; 

		// Update generation count
		generation_cnt++;

		#if defined (DEBUG_KRNL_GA)
		printf("eval_cnt: %u, generation_cnt: %u\n", eval_cnt, generation_cnt);
		#endif
	} // End while eval_cnt & generation_cnt

	// ------------------------------------------------------------------
	// Write final pop & energies back to FPGA-board DDRs
	// ------------------------------------------------------------------
	for (ushort pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
		for (uchar gene_cnt = 0; gene_cnt < DockConst_num_of_genes; gene_cnt++) {
			GlobPopCurrFinal[pop_cnt*ACTUAL_GENOTYPE_LENGTH + gene_cnt] = LocalPopCurr[pop_cnt][gene_cnt];
		}

		GlobEneCurr[pop_cnt] = LocalEneCurr[pop_cnt];
	}

	#if defined (DEBUG_KRNL_GA)
	printf("GA: %u %u\n", active, DockConst_pop_size -1);
	#endif

	#if defined (DEBUG_ACTIVE_KERNEL)
	printf("	%-20s: %s\n", "Krnl_GA", "disabled");
	#endif

	// Write final evals & generation counts to FPGA-board DDRs
	GlobEvals_performed[Host_RunId] = eval_cnt;
	GlobGens_performed [Host_RunId] = generation_cnt;

	/*printf("Host_RunId: %u, eval_cnt: %u, generation_cnt: %u\n", Host_RunId, eval_cnt, generation_cnt);*/

  	#ifdef PRINT_ALL_KRNL
  	printf("\tFinishing libkernel_ga\n");
  	#endif
	return 0;
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
