#include "defines.h"
#include "math.h"

//IC:  initial calculation of energy of populations
//GG:  genetic generation 
//LS:  local search
//OFF: turn off 

#define PIPE_DEPTH_16  16
#define PIPE_DEPTH_64  64
#define PIPE_DEPTH_512 512

// Status of pipe operation
// Success: 0
// Failure: negative value, e.g.: -1, -2, etc

// Important: the evaluation of failure of "pipe-expr" 
// must be done: (pipe-expr != PIPE_STATUS_SUCCESS),
// as a failure is characterize by any negative integer number.
typedef int nb_pipe_status;
#define PIPE_STATUS_SUCCESS      0
#define PIPE_STATUS_FAILURE	-1

// --------------------------------------------------------------------------
// Map the argument into the interval 0 - 180, or 0 - 360
// by adding/subtracting n*ang_max to/from it.
// Originally from: searchoptimum.c
// --------------------------------------------------------------------------
__attribute__((always_inline))
float map_angle_180(float angle)
{
	float x = angle;
	//while (x < 0.0f)
	if (x < 0.0f)   
	{ x += 180.0f; }
	//while (x > 180.0f)
	if (x > 180.0f) 
	{ x -= 180.0f; }
	return x;
}

__attribute__((always_inline))
float map_angle_360(float angle)
{
	float x = angle;
	//while (x < 0.0f)
	if (x < 0.0f)
	{ x += 360.0f; }
	//while (x > 360.0f)
	if (x > 360.0f)
	{ x -= 360.0f;}
	return x;
}

// --------------------------------------------------------------------------
// Lamarckian Genetic-Algorithm (GA): GA + LS (Local Search) 
// Originally from: searchoptimum.c
// --------------------------------------------------------------------------
void libkernel_ga (
	const 	float*           restrict GlobPopulationCurrentInitial,
			float*           restrict GlobPopulationCurrentFinal,
			float*           restrict GlobEnergyCurrent,
			unsigned int*    restrict GlobEvals_performed,
            unsigned int*    restrict GlobGens_performed,
			unsigned int              DockConst_pop_size,
		    unsigned int              DockConst_num_of_energy_evals,
			unsigned int              DockConst_num_of_generations,
		    float                     DockConst_tournament_rate,
			float                     DockConst_mutation_rate,
		    float                     DockConst_abs_max_dmov,
			float                     DockConst_abs_max_dang,
		    float                     Host_two_absmaxdmov,
			float                     Host_two_absmaxdang,
			float                     DockConst_crossover_rate,
			unsigned int              DockConst_num_of_lsentities,
			unsigned char             DockConst_num_of_genes,
	        unsigned short            Host_RunId,
			unsigned int 	      	  Host_Offset_Pop,
			unsigned int	      	  Host_Offset_Ene
)
{
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

	float LocalPopCurr[MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];
	float LocalEneCurr[MAX_POPSIZE];

	const float* GlobPopCurrInitial = & GlobPopulationCurrentInitial [Host_Offset_Pop];
	      float* GlobPopCurrFinal   = & GlobPopulationCurrentFinal   [Host_Offset_Pop];
	      float* GlobEneCurr        = & GlobEnergyCurrent     	     [Host_Offset_Ene];

	// ------------------------------------------------------------------
	// Initial Calculation (IC) of scores
	// ------------------------------------------------------------------
	// LOOP_FOR_GA_IC_OUTER
	for (unsigned short pop_cnt = 0; pop_cnt < DockConst_pop_size; pop_cnt++) {
		// Calculate energy
		const int tmp_int_zero = 0;
		write_pipe_block(pipe00ga2igl00ic00active, &tmp_int_zero);

		LOOP_FOR_GA_IC_INNER_WRITE_GENOTYPE:
		for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
			float tmp_ic;
			tmp_ic = GlobPopCurrInitial[pop_cnt*ACTUAL_GENOTYPE_LENGTH + gene_cnt];

			LocalPopCurr[pop_cnt][gene_cnt & MASK_GENOTYPE] = tmp_ic;
			write_pipe_block(pipe00ic2conf00genotype, &tmp_ic);	
		}

		#if defined (DEBUG_KRNL_IC)
		printf("\nIC - tx pop: %u", pop_cnt); 		
		#endif

		// Read energy
		float energyIA_IC_rx;
		float energyIE_IC_rx;

		nb_pipe_status intra_valid = PIPE_STATUS_FAILURE;
		nb_pipe_status inter_valid = PIPE_STATUS_FAILURE;	

		// LOOP_WHILE_GA_IC_INNER_READ_ENERGY
		while( (intra_valid != PIPE_STATUS_SUCCESS) || (inter_valid != PIPE_STATUS_SUCCESS)) {

			if (intra_valid != PIPE_STATUS_SUCCESS) {
				intra_valid = read_pipe(pipe00intrae2storeic00intrae, &energyIA_IC_rx);
			}
			else if (inter_valid != PIPE_STATUS_SUCCESS) {
				inter_valid = read_pipe(pipe00intere2storeic00intere, &energyIE_IC_rx);
			}
		}

		LocalEneCurr[pop_cnt] = energyIA_IC_rx + energyIE_IC_rx;

		#if defined (DEBUG_KRNL_IC)
		printf(", IC - rx pop: %u\n", pop_cnt); 		
		#endif
	}
	// ------------------------------------------------------------------

	unsigned int eval_cnt = DockConst_pop_size; // takes into account the IC evals

	unsigned int generation_cnt = 0;

	// LOOP_WHILE_GA_MAIN
	while ((eval_cnt < DockConst_num_of_energy_evals) && (generation_cnt < DockConst_num_of_generations)) {

		float LocalPopNext[MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];
		float LocalEneNext[MAX_POPSIZE];

		// ------------------------------------------------------------------
		// Genetic Generation (GG)
		// ------------------------------------------------------------------
		float loc_energies[MAX_POPSIZE];
		unsigned short best_entity = 0;

		// LOOP_FOR_GA_SHIFT
//		for (unsigned short pop_cnt=1; pop_cnt<DockConst_pop_size; pop_cnt++) {
		for (unsigned short pop_cnt=0; pop_cnt<DockConst_pop_size; pop_cnt++) {
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
		// LOOP_FOR_GA_OUTER_GLOBAL
		for (unsigned short new_pop_cnt = 1; new_pop_cnt < DockConst_pop_size; new_pop_cnt++) {

			// ---------------------------------------------------
			// Elitism: copying the best entity to new population
			// ---------------------------------------------------
			if (new_pop_cnt == 1) {
				// LOOP_FOR_GA_INNER_ELITISM
				for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
					LocalPopNext[0][gene_cnt & MASK_GENOTYPE] = LocalPopCurr[best_entity][gene_cnt & MASK_GENOTYPE]; 	
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

			// Get ushort binary_tournament selection prngs (parent index)
			// Get float binary_tournament selection prngs (tournament rate)
			float8 bt_tmp;
			read_pipe_block(pipe00prng2ga00bt00ushort00float00prng, &bt_tmp);

//printf("test point 1\n");
			// Convert: float prng that must be still converted to short
			float bt_tmp_uf0 = bt_tmp.s0;
			float bt_tmp_uf1 = bt_tmp.s2;
			float bt_tmp_uf2 = bt_tmp.s4;
			float bt_tmp_uf3 = bt_tmp.s6;

			// short prng ready to be used, replace ushort prng_BT_U[4];
/*
			ushort bt_tmp_u0 = *(uint*)&bt_tmp_uf0;
			ushort bt_tmp_u1 = *(uint*)&bt_tmp_uf1;
			ushort bt_tmp_u2 = *(uint*)&bt_tmp_uf2;
			ushort bt_tmp_u3 = *(uint*)&bt_tmp_uf3;
*/
			// Check "Krnl_Prng_BT_ushort_float"
			// To surpass error in hw_emu		
			unsigned short bt_tmp_u0 = bt_tmp_uf0;
			unsigned short bt_tmp_u1 = bt_tmp_uf1;
			unsigned short bt_tmp_u2 = bt_tmp_uf2;
			unsigned short bt_tmp_u3 = bt_tmp_uf3;

			// float prng ready to used, replace float prng_BT_F[4];
			float bt_tmp_f0 = bt_tmp.s1;
			float bt_tmp_f1 = bt_tmp.s3;
			float bt_tmp_f2 = bt_tmp.s5;
			float bt_tmp_f3 = bt_tmp.s7;

			unsigned short parent1;
			unsigned short parent2; 

			// First parent
			if (loc_energies[bt_tmp_u0] < loc_energies[bt_tmp_u1]) {
				if (bt_tmp_f0 < DockConst_tournament_rate) {parent1 = bt_tmp_u0;}
				else				           {parent1 = bt_tmp_u1;}}
			else {
				if (bt_tmp_f1 < DockConst_tournament_rate) {parent1 = bt_tmp_u1;}
				else				           {parent1 = bt_tmp_u0;}}

			// The better will be the second parent
			if (loc_energies[bt_tmp_u2] < loc_energies[bt_tmp_u3]) {
				if (bt_tmp_f2 < DockConst_tournament_rate) {parent2 = bt_tmp_u2;}
				else		          	           {parent2 = bt_tmp_u3;}}
			else {
				if (bt_tmp_f3 < DockConst_tournament_rate) {parent2 = bt_tmp_u3;}
				else			                   {parent2 = bt_tmp_u2;}}

			// LOOP_FOR_GA_INNER_BT
			// local_entity_1 and local_entity_2 are population-parent1, population-parent2
			for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
				local_entity_1[gene_cnt & MASK_GENOTYPE] = LocalPopCurr[parent1][gene_cnt & MASK_GENOTYPE];
				local_entity_2[gene_cnt & MASK_GENOTYPE] = LocalPopCurr[parent2][gene_cnt & MASK_GENOTYPE];
			}

			// ---------------------------------------------------
			// Mating parents
			// ---------------------------------------------------	

			// get uchar genetic_generation prngs (gene index)
			// get float genetic_generation prngs (mutation rate)
			uchar2 prng_GG_C;
			read_pipe_block(pipe00prng2ga00gg00uchar00prng, &prng_GG_C);

//printf("test point 2\n");

			unsigned char covr_point_low;
			unsigned char covr_point_high;
			bool twopoint_cross_yes = false;

			if (prng_GG_C.x == prng_GG_C.y) {covr_point_low = prng_GG_C.x;}
			else {
				twopoint_cross_yes = true;
				if (prng_GG_C.x < prng_GG_C.y) { covr_point_low  = prng_GG_C.x;
					                         covr_point_high = prng_GG_C.y; }
				else {		      		 covr_point_low  = prng_GG_C.y;
   								 covr_point_high = prng_GG_C.x; }
			}
			
			// Reuse of bt prng float as crossover-rate
			bool crossover_yes = (DockConst_crossover_rate > bt_tmp_f0);

			const int tmp_int_zero = 0;
			write_pipe_block(pipe00ga2igl00gg00active, &tmp_int_zero);

//printf("test point 3\n");

			// LOOP_FOR_GA_INNER_CROSS_MUT
			for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
				float prngGG;
				read_pipe_block(pipe00prng2ga00gg00float00prng, &prngGG);

//printf("test point 4\n");

				float tmp_offspring;

				// Performing crossover
				if (   	(
					crossover_yes && (										// crossover
					( (twopoint_cross_yes == true)  && ((gene_cnt <= covr_point_low) || (gene_cnt > covr_point_high)) )  ||	// two-point crossover 			 		
					( (twopoint_cross_yes == false) && (gene_cnt <= covr_point_low))  					// one-point crossover
					)) || 
					(!crossover_yes)	// no crossover
				   ) {
					tmp_offspring = local_entity_1[gene_cnt & MASK_GENOTYPE];
				}
				else {
					tmp_offspring = local_entity_2[gene_cnt & MASK_GENOTYPE];
				}

				// Performing mutation
				if (DockConst_mutation_rate > prngGG) {
					if(gene_cnt<3) {
						tmp_offspring = tmp_offspring + Host_two_absmaxdmov*prngGG-DockConst_abs_max_dmov;
					}
					else {
						float tmp;
						tmp = tmp_offspring + Host_two_absmaxdang*prngGG-DockConst_abs_max_dang;
						if (gene_cnt==4) { tmp_offspring = map_angle_180(tmp); }
						else             { tmp_offspring = map_angle_360(tmp); }
					}
				}

				// Calculate energy
				LocalPopNext [new_pop_cnt][gene_cnt & MASK_GENOTYPE] = tmp_offspring;
				write_pipe_block(pipe00gg2conf00genotype, &tmp_offspring);
//printf("test point 5\n");
			}

			#if defined (DEBUG_KRNL_GG)
			printf("GG - tx pop: %u", new_pop_cnt); 		
			#endif	

			// Read energy
			float energyIA_GG_rx;
			float energyIE_GG_rx;

			nb_pipe_status intra_valid = PIPE_STATUS_FAILURE;
			nb_pipe_status inter_valid = PIPE_STATUS_FAILURE;

			// LOOP_WHILE_GA_INNER_READ_ENERGIES
			while( (intra_valid != PIPE_STATUS_SUCCESS) || (inter_valid != PIPE_STATUS_SUCCESS)) {

				if (intra_valid != PIPE_STATUS_SUCCESS) {
					intra_valid = read_pipe(pipe00intrae2storegg00intrae, &energyIA_GG_rx);
				}
				else if (inter_valid != PIPE_STATUS_SUCCESS) {
					inter_valid = read_pipe(pipe00intere2storegg00intere, &energyIE_GG_rx);
				}

//printf("intra_valid: %i, inter_valid: %i\n", intra_valid, inter_valid);
			}
//printf("test point 5\n");			
			LocalEneNext[new_pop_cnt] = energyIA_GG_rx + energyIE_GG_rx;

			#if defined (DEBUG_KRNL_GG)
			printf(", GG - rx pop: %u\n", new_pop_cnt); 		
			#endif
		} 
		// ------------------------------------------------------------------
		// LS: Local Search
		// Subject num_of_entity_for_ls pieces of offsprings to LS 
		// ------------------------------------------------------------------

		unsigned int ls_eval_cnt = 0;

		/*
		#pragma ivdep
		*/
		// LOOP_FOR_GA_LS_OUTER
		for (unsigned short ls_ent_cnt=0; ls_ent_cnt<DockConst_num_of_lsentities; ls_ent_cnt+=9) {

			// Choose random & different entities on every iteration
			ushort16 entity_ls;
			read_pipe_block(pipe00prng2ga00ls12300ushort00prng, &entity_ls);

//printf("test point LS 1\n");

			unsigned short entity_ls1 = entity_ls.s0;
			unsigned short entity_ls2 = entity_ls.s1;
			unsigned short entity_ls3 = entity_ls.s2;
			unsigned short entity_ls4 = entity_ls.s3;
			unsigned short entity_ls5 = entity_ls.s4;
			unsigned short entity_ls6 = entity_ls.s5;
			unsigned short entity_ls7 = entity_ls.s6;
			unsigned short entity_ls8 = entity_ls.s7;
			unsigned short entity_ls9 = entity_ls.s8;

			write_pipe_block(pipe00ga2ls00ls100energy, &LocalEneNext[entity_ls1]);
			write_pipe_block(pipe00ga2ls00ls200energy, &LocalEneNext[entity_ls2]);
			write_pipe_block(pipe00ga2ls00ls300energy, &LocalEneNext[entity_ls3]);
			write_pipe_block(pipe00ga2ls00ls400energy, &LocalEneNext[entity_ls4]);
			write_pipe_block(pipe00ga2ls00ls500energy, &LocalEneNext[entity_ls5]);
			write_pipe_block(pipe00ga2ls00ls600energy, &LocalEneNext[entity_ls6]);
			write_pipe_block(pipe00ga2ls00ls700energy, &LocalEneNext[entity_ls7]);
			write_pipe_block(pipe00ga2ls00ls800energy, &LocalEneNext[entity_ls8]);
			write_pipe_block(pipe00ga2ls00ls900energy, &LocalEneNext[entity_ls9]);

//printf("test point LS 2\n");

			// LOOP_GA_LS_INNER_WRITE_GENOTYPE
			for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
				write_pipe_block(pipe00ga2ls00ls100genotype, &LocalPopNext[entity_ls1][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls200genotype, &LocalPopNext[entity_ls2][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls300genotype, &LocalPopNext[entity_ls3][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls400genotype, &LocalPopNext[entity_ls4][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls500genotype, &LocalPopNext[entity_ls5][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls600genotype, &LocalPopNext[entity_ls6][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls700genotype, &LocalPopNext[entity_ls7][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls800genotype, &LocalPopNext[entity_ls8][gene_cnt & MASK_GENOTYPE]);
				write_pipe_block(pipe00ga2ls00ls900genotype, &LocalPopNext[entity_ls9][gene_cnt & MASK_GENOTYPE]);
			}

//printf("test point LS 3\n");

			float2 evalenergy_tmp1;
			float2 evalenergy_tmp2;
			float2 evalenergy_tmp3;
			float2 evalenergy_tmp4;
			float2 evalenergy_tmp5;
			float2 evalenergy_tmp6;
			float2 evalenergy_tmp7;
			float2 evalenergy_tmp8;
			float2 evalenergy_tmp9;

			nb_pipe_status ls1_done = PIPE_STATUS_FAILURE;
			nb_pipe_status ls2_done = PIPE_STATUS_FAILURE;
			nb_pipe_status ls3_done = PIPE_STATUS_FAILURE;
		 	nb_pipe_status ls4_done = PIPE_STATUS_FAILURE;
			nb_pipe_status ls5_done = PIPE_STATUS_FAILURE;
			nb_pipe_status ls6_done = PIPE_STATUS_FAILURE;
			nb_pipe_status ls7_done = PIPE_STATUS_FAILURE;
			nb_pipe_status ls8_done = PIPE_STATUS_FAILURE;
			nb_pipe_status ls9_done = PIPE_STATUS_FAILURE;  

			// LOOP_WHILE_GA_LS_INNER_READ_ENERGIES
			while( (ls1_done != PIPE_STATUS_SUCCESS) || 
			       (ls2_done != PIPE_STATUS_SUCCESS) || 
			       (ls3_done != PIPE_STATUS_SUCCESS) || 
			       (ls4_done != PIPE_STATUS_SUCCESS) || 
			       (ls5_done != PIPE_STATUS_SUCCESS) ||
			       (ls6_done != PIPE_STATUS_SUCCESS) || 
			       (ls7_done != PIPE_STATUS_SUCCESS) || 
			       (ls8_done != PIPE_STATUS_SUCCESS) || 
			       (ls9_done != PIPE_STATUS_SUCCESS) 
			)
			{
				if (ls1_done != PIPE_STATUS_SUCCESS) {
					ls1_done = read_pipe(pipe00ls2ga00ls100evalenergy, &evalenergy_tmp1);
				}
				else if (ls2_done != PIPE_STATUS_SUCCESS) {
					ls2_done = read_pipe(pipe00ls2ga00ls200evalenergy, &evalenergy_tmp2);
				}
				else if (ls3_done != PIPE_STATUS_SUCCESS) {
					ls3_done = read_pipe(pipe00ls2ga00ls300evalenergy, &evalenergy_tmp3);
				}
				else if (ls4_done != PIPE_STATUS_SUCCESS) {
					ls4_done = read_pipe(pipe00ls2ga00ls400evalenergy, &evalenergy_tmp4);
				}
				else if (ls5_done != PIPE_STATUS_SUCCESS) {
					ls5_done = read_pipe(pipe00ls2ga00ls500evalenergy, &evalenergy_tmp5);
				}
				else if (ls6_done != PIPE_STATUS_SUCCESS) {
					ls6_done = read_pipe(pipe00ls2ga00ls600evalenergy, &evalenergy_tmp6);
				}
				else if (ls7_done != PIPE_STATUS_SUCCESS) {
					ls7_done = read_pipe(pipe00ls2ga00ls700evalenergy, &evalenergy_tmp7);
				}
				else if (ls8_done != PIPE_STATUS_SUCCESS) {
					ls8_done = read_pipe(pipe00ls2ga00ls800evalenergy, &evalenergy_tmp8);
				}
				else if (ls9_done != PIPE_STATUS_SUCCESS) {
					ls9_done = read_pipe(pipe00ls2ga00ls900evalenergy, &evalenergy_tmp9);
				}
			}
		
			#if defined (DEBUG_KRNL_LS)
			printf("LS - got all eval & energies back\n");
			#endif

			float eetmp1 = evalenergy_tmp1.x;
			float eetmp2 = evalenergy_tmp2.x;
			float eetmp3 = evalenergy_tmp3.x;
			float eetmp4 = evalenergy_tmp4.x;
			float eetmp5 = evalenergy_tmp5.x;
			float eetmp6 = evalenergy_tmp6.x;
			float eetmp7 = evalenergy_tmp7.x;
			float eetmp8 = evalenergy_tmp8.x;
			float eetmp9 = evalenergy_tmp9.x;

			unsigned int eval_tmp1 = *(unsigned int*)&eetmp1;
			unsigned int eval_tmp2 = *(unsigned int*)&eetmp2;
			unsigned int eval_tmp3 = *(unsigned int*)&eetmp3;
			unsigned int eval_tmp4 = *(unsigned int*)&eetmp4;
			unsigned int eval_tmp5 = *(unsigned int*)&eetmp5;
			unsigned int eval_tmp6 = *(unsigned int*)&eetmp6;
			unsigned int eval_tmp7 = *(unsigned int*)&eetmp7;
			unsigned int eval_tmp8 = *(unsigned int*)&eetmp8;
			unsigned int eval_tmp9 = *(unsigned int*)&eetmp9;

			LocalEneNext[entity_ls1] = evalenergy_tmp1.y;
			LocalEneNext[entity_ls2] = evalenergy_tmp2.y;
			LocalEneNext[entity_ls3] = evalenergy_tmp3.y;
			LocalEneNext[entity_ls4] = evalenergy_tmp4.y;
			LocalEneNext[entity_ls5] = evalenergy_tmp5.y;
			LocalEneNext[entity_ls6] = evalenergy_tmp6.y;
			LocalEneNext[entity_ls7] = evalenergy_tmp7.y;
			LocalEneNext[entity_ls8] = evalenergy_tmp8.y;
			LocalEneNext[entity_ls9] = evalenergy_tmp9.y;

			/*
			#pragma ivdep
			*/
			// LOOP_FOR_GA_LS_INNER_READ_GENOTYPE
			for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {

				read_pipe_block(pipe00ls2ga00ls100genotype, &LocalPopNext[entity_ls1][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls200genotype, &LocalPopNext[entity_ls2][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls300genotype, &LocalPopNext[entity_ls3][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls400genotype, &LocalPopNext[entity_ls4][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls500genotype, &LocalPopNext[entity_ls5][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls600genotype, &LocalPopNext[entity_ls6][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls700genotype, &LocalPopNext[entity_ls7][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls800genotype, &LocalPopNext[entity_ls8][gene_cnt & MASK_GENOTYPE]);
				read_pipe_block(pipe00ls2ga00ls900genotype, &LocalPopNext[entity_ls9][gene_cnt & MASK_GENOTYPE]);
			}

			ls_eval_cnt += eval_tmp1 + eval_tmp2 + eval_tmp3 + eval_tmp4 + eval_tmp5 + eval_tmp6 + eval_tmp7 + eval_tmp8 + eval_tmp9;

			#if defined (DEBUG_KRNL_LS)
			printf("%u, ls_eval_cnt: %u\n", ls_ent_cnt, ls_eval_cnt);
			printf("LS - got all genotypes back\n");
			#endif
		} // End of for-loop ls_ent_cnt
		// ------------------------------------------------------------------

		// Update current pops & energies
		// LOOP_FOR_GA_UPDATEPOP_OUTER
		for (unsigned short pop_cnt=0; pop_cnt<DockConst_pop_size; pop_cnt++) {
			// LOOP_GA_UPDATEPOP_INNER
			for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
				LocalPopCurr[pop_cnt][gene_cnt & MASK_GENOTYPE] = LocalPopNext[pop_cnt][gene_cnt & MASK_GENOTYPE];
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
	// Off: turn off all other kernels
	// ------------------------------------------------------------------

	// Turn off PRNG kernels
	const int tmp_int_one = 1;
	write_pipe_block(pipe00ga2prng00bt00ushort00float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00gg00uchar00off, 	&tmp_int_one);
	write_pipe_block(pipe00ga2prng00gg00float00off, 	&tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls12300ushort00off,  	&tmp_int_one);

	write_pipe_block(pipe00ga2prng00ls00float00off,  &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls200float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls300float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls400float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls500float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls600float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls700float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls800float00off, &tmp_int_one);
	write_pipe_block(pipe00ga2prng00ls900float00off, &tmp_int_one);

	// Turn off LS kernels
	write_pipe_block(pipe00ga2ls00off100active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off200active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off300active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off400active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off500active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off600active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off700active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off800active, &tmp_int_one);
	write_pipe_block(pipe00ga2ls00off900active, &tmp_int_one);

	// Turn off IGL_Arbiter, Conform, InterE, IntraE kernerls
	write_pipe_block(pipe00iglarbiter00off,     		&tmp_int_one);

	// Write final pop & energies back to FPGA-board DDRs
	// LOOP_GA_WRITEPOP2DDR_OUTER
	for (unsigned short pop_cnt=0;pop_cnt<DockConst_pop_size; pop_cnt++) { 	
		// LOOP_GA_WRITEPOP2DDR_INNER
		for (unsigned char gene_cnt=0; gene_cnt<DockConst_num_of_genes; gene_cnt++) {
			GlobPopCurrFinal[pop_cnt*ACTUAL_GENOTYPE_LENGTH + gene_cnt] = LocalPopCurr[pop_cnt][gene_cnt & MASK_GENOTYPE];
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
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

/*
#include "Krnl_PRNG.cl"

#include "Krnl_LS.cl"
#include "Krnl_LS2.cl"
#include "Krnl_LS3.cl"
#include "Krnl_LS4.cl"
#include "Krnl_LS5.cl"
#include "Krnl_LS6.cl"
#include "Krnl_LS7.cl"
#include "Krnl_LS8.cl"
#include "Krnl_LS9.cl"

#include "Krnl_IGL_SimplifArbiter.cl"

#include "Krnl_Conform.cl"
#include "Krnl_InterE.cl"
#include "Krnl_IntraE.cl"
*/
