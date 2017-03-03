// Enable the channels extension
#pragma OPENCL EXTENSION cl_altera_channels : enable

// Define kernel file-scope channel variable
// Buffered channels ACTUAL_GENOTYPE_LENGTH	(MAX_NUM_OF_ROTBONDS+6) =38
channel float chan_GA2Conf_genotype __attribute__((depth(38)));

// Buffered channels MAX_NUM_OF_ATOMS=90
//channel float chan_GA2Conf_ligandatom_idxyzq		__attribute__((depth(1280)));
channel float chan_GA2Conf_x	__attribute__((depth(90)));
channel float chan_GA2Conf_y	__attribute__((depth(90)));
channel float chan_GA2Conf_z	__attribute__((depth(90)));

//channel float chan_Conf2Intere_ligandatom_idxyzq	__attribute__((depth(1280)));
//channel float chan_Conf2Intrae_ligandatom_idxyzq	__attribute__((depth(1280)));
channel float chan_Conf2Intere_x __attribute__((depth(90)));
channel float chan_Conf2Intere_y __attribute__((depth(90)));
channel float chan_Conf2Intere_z __attribute__((depth(90)));

channel float chan_Conf2Intrae_x __attribute__((depth(90)));
channel float chan_Conf2Intrae_y __attribute__((depth(90)));
channel float chan_Conf2Intrae_z __attribute__((depth(90)));

channel float chan_Intere2GA_intere;
channel float chan_Intrae2GA_intrae;

#include "auxiliary_genetic.cl"
#include "auxiliary_performls.cl"

// Next structures were copied from calcenergy.h
typedef struct
{
	char  	 	num_of_atoms;
	char   		num_of_atypes;
	int    		num_of_intraE_contributors;
	char   		gridsize_x;
	char   		gridsize_y;
	char   		gridsize_z;
	float  		grid_spacing;
/*
	float* 		fgrids;
*/
	int    		rotbondlist_length;
	float  		coeff_elec;
	float  		coeff_desolv;
/*
	float* 		conformations_current;
	float* 		energies_current;
	float* 		conformations_next;
	float* 		energies_next;
	int*   		evals_of_new_entities;
	unsigned int* 	prng_states;
*/
	int    		pop_size;
	int    		num_of_genes;
	float  		tournament_rate;
	float  		crossover_rate;
	float  		mutation_rate;
	float  		abs_max_dmov;
	float  		abs_max_dang;
	float  		lsearch_rate;
	unsigned int 	num_of_lsentities;
	float  		rho_lower_bound;
	float  		base_dmov_mul_sqrt3;
	float  		base_dang_mul_sqrt3;
	unsigned int 	cons_limit;
	unsigned int 	max_num_of_iters;
	float  		qasp;
} Dockparameters;

// Constant struct
typedef struct
{
       float atom_charges_const[MAX_NUM_OF_ATOMS];
       char  atom_types_const  [MAX_NUM_OF_ATOMS];
       char  intraE_contributors_const[3*MAX_INTRAE_CONTRIBUTORS];
       float VWpars_AC_const   [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
       float VWpars_BD_const   [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
       float dspars_S_const    [MAX_NUM_OF_ATYPES];
       float dspars_V_const    [MAX_NUM_OF_ATYPES];
       int   rotlist_const     [MAX_NUM_OF_ROTATIONS];
       float ref_coords_x_const[MAX_NUM_OF_ATOMS];
       float ref_coords_y_const[MAX_NUM_OF_ATOMS];
       float ref_coords_z_const[MAX_NUM_OF_ATOMS];
       float rotbonds_moving_vectors_const[3*MAX_NUM_OF_ROTBONDS];
       float rotbonds_unit_vectors_const  [3*MAX_NUM_OF_ROTBONDS];
       float ref_orientation_quats_const  [4*MAX_NUM_OF_RUNS];
} kernelconstant;





// --------------------------------------------------------------------------
// The function performs a generational genetic algorithm based search 
// on the search space.
// The first parameter is the population which must be filled with initial values 
// before calling this function. 
// The other parameters are variables which describe the grids, 
// the docking parameters and the ligand to be docked. 
// Originally from: searchoptimum.c
// --------------------------------------------------------------------------
__kernel __attribute__ ((reqd_work_group_size(1,1,1)))
/*
void Krnl_GA(__global       float*           restrict GlobPopulation,
	     __global       float*           restrict GlobLigandAtom_idxyzq,
             __global const float*           restrict GlobFgrids,
	     __global       uint*            restrict GlobCounter,
	     __global       uint*            restrict GlobPRNG,
	     __global const Dockingconstant* restrict DockConst,
	     __global const Ligandconstant*  restrict LigConst,
	     __global const Gridconstant*    restrict GridConst)
*/
void Krnl_GA(__global const float*           restrict GlobFgrids,
	     __global       float*           restrict GlobPopulationCurrent,
	     __global       float*           restrict GlobEnergyCurrent,
	     __global       float*           restrict GlobPopulationNext,
	     __global       float*           restrict GlobEnergyNext,
             __global const float*           restrict GlobPRNG,
	     __global const kernelconstant* restrict KerConst,
	     __global const Dockparameters* restrict DockConst)

{
#ifdef EMULATOR
	printf("Krnl GA!!\n");
#endif

	uint parent1, parent2;
	uint entity_for_ls;
	uint i;
	uint j;
	uint LS_eval;
	uint new_pop_cnt;
	uint best_entity_id;
	uint i_pop;

	// --------------------------------------
	// L30nardoSV
	__local float local_entity[40];		//float private_entity[40];
	__local float local_entity_2[40];	//float private_entity_2[40];
	
	// The second is required only because gen_new_genotype() returns two offsprings
	__local float offspring1_genotype [40];

	// Local mem to store new_population
	__local float loc_new_population [MAX_POPSIZE*40];

	//Liganddata myligand_temp;
	__local float loc_templigand_atom_idxyzq [MAX_NUM_OF_ATOMS*5];
	// --------------------------------------

	float avg_energy;
	uint evals_for_ls_in_this_cycle;

	uint eval_cnt = 0;
	uint generation_cnt = 1;

	uint num_of_evals_for_ls  = 0;

	//capturing and converting algorithm parameters

	#if defined (DEBUG_KERNEL1)
	printf("\nParameters of the genetic algorihtm:\n");
	printf("Rate of crossover: %f%%\n",crossover_rate);
	printf("Rate of mutation: %f%%\n",mutation_rate);
	printf("Rate of local search: %f%%\n",lsearch_rate);
	printf("Population size: %d\n", pop_size);
	printf("Maximal delta movement during mutation: +/-%fA\n",abs_max_dmov);
	printf("maximal delta angle during mutation: +/-%f°\n", abs_max_dang);
	printf("Rho lower bound: %f\n",rho_lower_bound); 
	printf("Maximal delta movement during ls: +/-%fA\n",base_dmov_mul_sqrt3);
	printf("Maximal delta angle during ls: +/-%f°\n",base_dang_mul_sqrt3);
	#endif

	#if defined (DEBUG_KERNEL1)
	printf("\nmyligand_num_of_atoms: %u\n", myligand_num_of_atoms);
	printf("myligand_num_of_rotbonds: %u\n",  myligand_num_of_rotbonds);
	printf("\nmyginfo_size_x: %u\n", myginfo_size_x);
	printf("myginfo_size_y: %u\n",   myginfo_size_y);
	printf("myginfo_size_z: %u\n",   myginfo_size_z);
	printf("myginfo_spacing: %f\n",  myginfo_spacing);
	printf("myginfo_num_of_atypes: %u\n",  myginfo_num_of_atypes);
	printf("\nmypars_num_of_energy_evals: %u\n",  mypars_num_of_energy_evals);
	printf("mypars_num_of_generations: %u\n",  mypars_num_of_generations);
	printf("mypars_dmov_mask: %u\n",  mypars_dmov_mask);
	printf("mypars_dang_mask: %u\n",  mypars_dang_mask);
	printf("mypars_mutation_rate: %u\n",  mypars_mutation_rate);
	printf("mypars_crossover_rate: %u\n",  mypars_crossover_rate);
	printf("mypars_lsearch_rate: %f\n",  mypars_lsearch_rate);
	printf("mypars_tournament_rate: %u\n",  mypars_tournament_rate);
	printf("mypars_rho_lower_bound: %u\n",  mypars_rho_lower_bound);
	printf("mypars_base_dmov_mul_sqrt3: %u\n",  mypars_base_dmov_mul_sqrt3);
	printf("mypars_base_dang_mul_sqrt3: %u\n",  mypars_base_dang_mul_sqrt3);
	printf("mypars_cons_limit: %u\n",  mypars_cons_limit);
	printf("mypars_max_num_of_iters: %u\n",  mypars_max_num_of_iters);
	printf("mypars_pop_size: %u\n",  mypars_pop_size);
	#endif

	//calculating energies of initial population

	uint pop_cnt;



	LOOP_GEN_GENERATIONAL_1:
	for (pop_cnt = 0; pop_cnt < DockConst->pop_size; pop_cnt++)
	{	
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Krnl_Conform
		for (uint pipe_cnt=0; pipe_cnt<DockConst->num_of_atoms; pipe_cnt++) {
			write_channel_altera(chan_GA2Conf_x, KerConst->ref_coords_x_const[pipe_cnt]);
			write_channel_altera(chan_GA2Conf_y, KerConst->ref_coords_y_const[pipe_cnt]);
			write_channel_altera(chan_GA2Conf_z, KerConst->ref_coords_z_const[pipe_cnt]);
		}


		for (uint pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
			write_channel_altera(chan_GA2Conf_genotype, GlobPopulationCurrent[pop_cnt*ACTUAL_GENOTYPE_LENGTH+ pipe_cnt]);
		}


		////Krnl_InterE
		//GlobPopulation[init_cnt*40+39] = read_channel_altera(chan_Intere2GA_intere);
		
		////Krnl_IntraE
		//GlobPopulation[init_cnt*40+38] = read_channel_altera(chan_Intrae2GA_intrae);

		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
		#if defined (DEBUG_LEO)
		printf("pop_cnt (INC): % u\n", pop_cnt);
		#endif
	} // End of LOOP_GEN_GENERATIONAL_1
	
	eval_cnt = pop_cnt++;
	
	#if defined (DEBUG_LEO)
	printf("eval_cnt (INC): % u\n", eval_cnt);
	#endif




	// REQUIRES FENCES TO SYNC 

/*



	// find_best
	// to store energies read from populations
	__local float loc_pop_i[MAX_POPSIZE*2];
	__local float loc_pop_energies[MAX_POPSIZE];

	// binary_tournament
	__local float loc_pop_par0[2];
	__local float loc_pop_par1[2];

	// local search
	__local float entity_possible_new_genotype [40];
	__local float loc_LS_templigand_atom_idxyzq [MAX_NUM_OF_ATOMS*5];




	// **********************************************
	// ADD VENDOR SPECIFIC PRAGMA
	// **********************************************
	LOOP_WHILE_GEN_GENERATIONAL_1:
	while ((eval_cnt < DockConst->num_of_energy_evals) && (generation_cnt < DockConst->num_of_generations))
	{
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


		//creating a new population

		//identifying best entity
		best_entity_id = find_best(GlobPopulation, 
					   loc_pop_i,
					   loc_pop_energies,
					   DockConst->pop_size); 

		#if defined (DEBUG_KERNEL1)
		avg_energy = 0.0f;
		printf("\n\n\nFinal state of the %u. generation:\n", generation_cnt);
		printf("----------------------------\n\n");
		for (i=0; i<DockConst->pop_size; i++)
		{
			avg_energy += GlobPopulation [i*40 + 38] + GlobPopulation [i*40 + 39];
			printf("Entity %3u: ", i);
			for (j=0; j<myligand_num_of_rotbonds+6; j++) {
				printf("%8.3f ", GlobPopulation [i*40 +j]);
			}
			printf("   energies: %10.3f %10.3f (sum: %10.3f)\n", 
				   GlobPopulation [i*40 + 38],
				   GlobPopulation [i*40 + 39], 
				   GlobPopulation [i*40 + 38] + GlobPopulation [i*40 + 39]);
		}
		printf("\nAverage energy: %f\nBest energy: %f %f (sum: %f)\n\n", avg_energy/DockConst->pop_size, 
		       GlobPopulation[best_entity_id*40 + 38], 
		       GlobPopulation[best_entity_id*40 + 39],
		       GlobPopulation[best_entity_id*40 + 38] + GlobPopulation[best_entity_id*40 + 39]);
		#endif

		//elitism - copying the best entity to new population
		async_work_group_copy(loc_new_population, GlobPopulation+best_entity_id*40, 40, 0);


		//new population consists of one member currently
		//new_pop_cnt = 1;









//







		// -----------------------------------------------------------------------
		// **********************************************
		// ADD VENDOR SPECIFIC PRAGMA
		// **********************************************
		LOOP_GEN_GENERATIONAL_2:
		for (new_pop_cnt = 1; new_pop_cnt < DockConst->pop_size; new_pop_cnt++)
		{
			//selecting two individuals randomly
			binary_tournament_selection(GlobPopulation,
						    DockConst->pop_size,
						    &parent1,
						    &parent2,
						    DockConst->tournament_rate,
						    GlobPRNG,
						    loc_pop_par0,
						    loc_pop_par1);

			//mating parents	
			async_work_group_copy(local_entity,   GlobPopulation+parent1*40, 40, 0);
			async_work_group_copy(local_entity_2, GlobPopulation+parent2*40, 40, 0);

			// first two args are population [parent1], population [parent2]
			gen_new_genotype(GlobPRNG,
                                         local_entity, 
					 local_entity_2,
					 DockConst->mutation_rate,
					 DockConst->dmov_mask,
					 DockConst->dang_mask,
					 DockConst->crossover_rate,
					 LigConst->num_of_rotbonds+6,
					 offspring1_genotype);

			//evaluating first offspring
			async_work_group_copy(loc_templigand_atom_idxyzq, GlobLigandAtom_idxyzq, MAX_NUM_OF_ATOMS*5, 0);

			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// Krnl_Conform
			uint pipe_cnt;

			for (pipe_cnt=0; pipe_cnt<LigConst->num_of_atoms*5; pipe_cnt++)
			{
				write_channel_altera(chan_GA2Conf_ligandatom_idxyzq, loc_templigand_atom_idxyzq[pipe_cnt]);
			}

			for (pipe_cnt=0; pipe_cnt<40; pipe_cnt++)
			{
				write_channel_altera(chan_GA2Conf_genotype, offspring1_genotype[pipe_cnt]);
			}


			//Krnl_InterE
			offspring1_genotype[39] = read_channel_altera(chan_Intere2GA_intere);
		
			//Krnl_IntraE
			offspring1_genotype[38] = read_channel_altera(chan_Intrae2GA_intrae);


			//copying first offspring to population
			// **********************************************
			// ADD VENDOR SPECIFIC PRAGMA
			// **********************************************
			LOOP_GEN_GENERATIONAL_3:
			for (i_pop = 0; i_pop < 38; i_pop++) {
				loc_new_population[new_pop_cnt*40 + i_pop] = offspring1_genotype[i_pop];
			}
	
			#if defined (DEBUG_LEO)
			printf("eval_cnt (INC): %u, new_pop_cnt (INC): %u\n", eval_cnt, new_pop_cnt);
			#endif
					
		} // End of LOOP_GEN_GENERATIONAL_2
		// -----------------------------------------------------------------------













//





























//












		#if defined (DEBUG_LEO)
		printf("End of loop of new_pop_cnt, new_pop_cnt = %u\n", new_pop_cnt);
		#endif

		//updating old population with new one
		// the last two values are not updated here as
		// these are updated in the intere and intrae kernels
		// CAN be merged with the previos loop
		// REQUIRES FENCES TO SYNC 
		async_work_group_copy(GlobPopulation, loc_new_population, MAX_POPSIZE*40, 0);



		#if defined (DEBUG_KERNEL1)
		best_entity_id = find_best(GlobPopulation, DockConst->pop_size); 
		find_best(DockConst->pop_size); 
		avg_energy = 0;
		printf("\n\n\nState of the %u. generation before local search:\n", generation_cnt+1);
		printf("----------------------------\n\n");
		for (i=0; i<DockConst->pop_size; i++)
		{
			avg_energy += GlobPopulation [i*40 + 38] + GlobPopulation [i*40 + 39];
			printf("Entity %3u: ", i);
			for (j=0; j<myligand_num_of_rotbonds+6; j++) {
				printf("%f ", GlobPopulation [i*40 + j]);
			}
			printf("   energies: %f %f (sum: %f)\n",
				   GlobPopulation [i*40 + 38], 
				   GlobPopulation [i*40 + 39],
				   GlobPopulation [i*40 + 38] + GlobPopulation [i*40 + 39]);
		}
		printf("\nAverage energy: %f\nBest energy: %f %f (sum: %f)\n\n", avg_energy/DockConst->pop_size,
			GlobPopulation [best_entity_id*40 + 38],
			GlobPopulation [best_entity_id*40 + 39],
			GlobPopulation [best_entity_id*40 + 38] + GlobPopulation [best_entity_id*40 + 39]);
		#endif

		evals_for_ls_in_this_cycle = 0;

		#if defined (DEBUG_LEO)
		printf("num_of_entity_for_ls: %u \n", DockConst->num_of_entity_for_ls);
		#endif






















		// -----------------------------------------------------------------------
		// **********************************************
		// ADD VENDOR SPECIFIC PRAGMA
		// **********************************************
		LOOP_GEN_GENERATIONAL_4:
		//subjecting num_of_entity_for_ls pieces of offsprings to LS
		for (i=0; i<DockConst->num_of_entity_for_ls; i++)	
		{
			//choosing an entity randomly,
			//and without checking if it has already been subjected to LS in this cycle
			entity_for_ls = myrand_uint(GlobPRNG, DockConst->pop_size);

			#if defined (DEBUG_LEO)
			printf("entity_for_ls: %u\n", entity_for_ls);
			#endif

			#if defined (DEBUG_KERNEL1)
			printf("Entity %u before local search: ", entity_for_ls);
			for (j=0; j<myligand_num_of_rotbonds+6; j++) {
				printf("%f ", GlobPopulation [entity_for_ls*40 + j]);
			}
			printf("   energies: %f %f\n", 
				   GlobPopulation [entity_for_ls*40 + 38],
				   GlobPopulation [entity_for_ls*40 + 39]);
			#endif

			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//PIPELINE INSIDE

			//performing local search
			async_work_group_copy(local_entity, GlobPopulation+entity_for_ls*40, 40, 0);

			perform_LS(GlobLigandAtom_idxyzq,
				   local_entity,
				   entity_possible_new_genotype,
				   loc_LS_templigand_atom_idxyzq,
				   // ligand params
				   LigConst->num_of_atoms,
				   LigConst->num_of_rotbonds,
				   LigConst->num_of_intraE_contributors, 
				   // grid params
				   GlobFgrids,
				   GridConst->size_x,
				   GridConst->size_y,
				   GridConst->size_z,
				   GridConst->g1,
				   GridConst->g2,
				   GridConst->g3,
				   GridConst->spacing,
				   GridConst->num_of_atypes, 
				   // dockpars params
				   DockConst->rho_lower_bound, 
				   DockConst->base_dmov_mul_sqrt3,
				   DockConst->base_dang_mul_sqrt3, 
				   DockConst->max_num_of_iters,
				   DockConst->cons_limit,
				   DockConst->cons_limit,  
				   &LS_eval,
				   GlobPRNG
				   );

			async_work_group_copy(GlobPopulation+entity_for_ls*40, local_entity, 40, 0);
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			eval_cnt += LS_eval;
			evals_for_ls_in_this_cycle += LS_eval;


			#if defined (DEBUG_KERNEL1)
			printf("Entity %u after local search (%u evaluations): ", entity_for_ls, LS_eval);
			for (j=0; j<myligand_num_of_rotbonds+6; j++) {printf("%f ", GlobPopulation [entity_for_ls*40 + j]);}
			printf("   energies: %f %f\n",
				   GlobPopulation [entity_for_ls*40 + 38],
				   GlobPopulation [entity_for_ls*40 + 39]);
			#endif

			num_of_evals_for_ls += LS_eval;

			#if defined (DEBUG_LEO)
			printf("eval_cnt (INC): %u, evals_for_ls_in_this_cycle(INC): % u\n", eval_cnt, evals_for_ls_in_this_cycle);
			#endif

		} // End of LOOP_GEN_GENERATIONAL_4

		// -----------------------------------------------------------------------

		generation_cnt++;

		#if defined (DEBUG_LEO)
		printf("generation_cnt (INC): %u\n", generation_cnt);
		#endif










//










	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	} // End of while ((eval_cnt < mypars->num_of_energy_evals) && (generation_cnt < mypars->num_of_generations))




*/





/*
	GlobCounter[0] = eval_cnt;
	GlobCounter[1] = generation_cnt;
*/


	#if defined (DEBUG_LEO)
	printf("eval_cnt: %u\n",  GlobCounter[0]);
	printf("generations_cnt: %u\n",  GlobCounter[1]);
	#endif

	#if defined (DEBUG_KERNEL1)
	printf("Energy evaluations for LS: %u out of %u\n", num_of_evals_for_ls, eval_cnt);
	printf("Number of generations: %u\n", generation_cnt);
	#endif
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------


#include "Krnl_Conform.cl"
#include "Krnl_InterE.cl"
#include "Krnl_IntraE.cl"


