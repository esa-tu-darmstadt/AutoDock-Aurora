// Enable the channels extension
#pragma OPENCL EXTENSION cl_altera_channels : enable

// Define kernel file-scope channel variable
// Buffered channels 
// MAX_NUM_OF_ATOMS=90
// ACTUAL_GENOTYPE_LENGTH (MAX_NUM_OF_ROTBONDS+6) =38
channel float chan_GA2Conf_genotype __attribute__((depth(38)));
//channel char  chan_GA2Manager_mode;	// mode 1 or I: init calculation energy, 2 or G: genetic generation, 3 or P: local search
//channel uint  chan_GA2Manager_cnt;	// population count
channel char  chan_GA2Conf_active;	// active 1: receiving Kernel is active, 0 receiving Kernel is disabled
channel char  chan_GA2Conf_mode;	// mode 1 or I: init calculation energy, 2 or G: genetic generation, 3 or P: local search
channel uint  chan_GA2Conf_cnt;		// population count

channel float chan_Conf2Intere_x __attribute__((depth(90)));
channel float chan_Conf2Intere_y __attribute__((depth(90)));
channel float chan_Conf2Intere_z __attribute__((depth(90)));
channel char  chan_Conf2Intere_active;	
channel char  chan_Conf2Intere_mode;	
channel uint  chan_Conf2Intere_cnt;	

channel float chan_Conf2Intrae_x __attribute__((depth(90)));
channel float chan_Conf2Intrae_y __attribute__((depth(90)));
channel float chan_Conf2Intrae_z __attribute__((depth(90)));
channel char  chan_Conf2Intrae_active;	
channel char  chan_Conf2Intrae_mode;	
channel uint  chan_Conf2Intrae_cnt;

channel float chan_Intere2Store_intere;
channel char  chan_Intere2Store_active;	
channel char  chan_Intere2Store_mode;	
channel uint  chan_Intere2Store_cnt;

channel float chan_Intrae2Store_intrae;
channel char  chan_Intrae2Store_active;	
channel char  chan_Intrae2Store_mode;	
channel uint  chan_Intrae2Store_cnt;

channel char  chan_Store2GA_ack;
channel float chan_Store2GA_LSenergy;

#include "../defines.h"

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

	// L30nardoSV added
	unsigned int num_of_energy_evals;
	unsigned int num_of_generations;

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
       //float ref_orientation_quats_const  [4*MAX_NUM_OF_RUNS];
       float ref_orientation_quats_const  [4];
} kernelconstant;


#include "auxiliary_genetic.cl"
#include "auxiliary_performls.cl"


// --------------------------------------------------------------------------
// The function performs a generational genetic algorithm based search 
// on the search space.
// The first parameter is the population which must be filled with initial values 
// before calling this function. 
// The other parameters are variables which describe the grids, 
// the docking parameters and the ligand to be docked. 
// Originally from: searchoptimum.c
// --------------------------------------------------------------------------
__kernel
void Krnl_GA(//__global const float*           restrict GlobFgrids,
	     __global       float*           restrict GlobPopulationCurrent,
	     __global       float*           restrict GlobEnergyCurrent,
	     __global       float*           restrict GlobPopulationNext,
	     __global       float*           restrict GlobEnergyNext,
             __global       unsigned int*    restrict GlobPRNG,
	     __global const kernelconstant*  restrict KerConst,
	     __global const Dockparameters*  restrict DockConst,
	     __global       unsigned int*    restrict GlobEvals_performed,
	     __global       unsigned int*    restrict GlobGenerations_performed)
{
	//Print algorithm parameters

	#if defined (DEBUG_KERNEL1)
	printf("\nParameters of the genetic algorihtm:\n");
	printf("\nLigand num_of_atoms: %u\n",  DockConst->num_of_atoms);
	printf("Ligand num_of_atypes:  %u\n",  DockConst->num_of_atypes);
	printf("Ligand num_of_intraE_contributors: %u\n",  DockConst->num_of_atypes);
	printf("Grid size_x: %u\n", 		DockConst->gridsize_x);
	printf("Grid size_y: %u\n",   		DockConst->gridsize_y);
	printf("Grid size_z: %u\n",   		DockConst->gridsize_z);
	printf("Grid spacing: %f\n",  		DockConst->grid_spacing);
	printf("Ligand rotbondlist_length: %u\n",  DockConst->rotbondlist_length);
	printf("Ligand coeff_elec: %f\n",  	DockConst->coeff_elec);
	printf("Ligand coeff_desolv: %f\n",  	DockConst->coeff_desolv);
	printf("\nnum_of_energy_evals: %u\n",   DockConst->num_of_energy_evals);
	printf("num_of_generations: %u\n",   	DockConst->num_of_generations);
	printf("Population size: %u\n",         DockConst->pop_size);
	printf("Number of genes: %u\n",         DockConst->num_of_genes);
	printf("Tournament rate: %f\n",  	DockConst->tournament_rate);
	printf("Crossover rate: %f\n",  	DockConst->crossover_rate);
	printf("Mutation rate: %f\n",  		DockConst->mutation_rate);
	printf("Maximal delta movement during mutation: +/-%fA\n", DockConst->abs_max_dmov);
	printf("maximal delta angle during mutation: +/-%f°\n",    DockConst->abs_max_dang);
	printf("LS rate: %f\n",  		DockConst->lsearch_rate);
	printf("LS num_of_lsentities: %u\n",    DockConst->num_of_lsentities);
	printf("LS rho_lower_bound: %f\n",      DockConst->rho_lower_bound);	 //Rho lower bound
	printf("LS base_dmov_mul_sqrt3: %f\n",  DockConst->base_dmov_mul_sqrt3); //Maximal delta movement during ls
	printf("LS base_dang_mul_sqrt3: %f\n",  DockConst->base_dang_mul_sqrt3); //Maximal delta angle during ls
	printf("LS cons_limit: %u\n",           DockConst->cons_limit);
	printf("LS max_num_of_iters: %u\n",     DockConst->max_num_of_iters);
	printf("qasp: %f\n",     DockConst->qasp);
	#endif







	//Calculating energies of initial population
	uint eval_cnt = 0;
	uint generation_cnt = 1;
	uint pop_cnt;
	uint new_pop_cnt;

	char active = 0;
	char mode   = 0;
	char ack    = 0;








///*









	LOOP_GEN_GENERATIONAL_1:
	for (pop_cnt = 0; pop_cnt < DockConst->pop_size; pop_cnt++)
	{	
		// --------------------------------------------------------------
		// Send genotypes to channel
		// --------------------------------------------------------------
		for (uint pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
			write_channel_altera(chan_GA2Conf_genotype, GlobPopulationCurrent[pop_cnt*ACTUAL_GENOTYPE_LENGTH+ pipe_cnt]);
		}
		active = 1;
		mode   = 1;
		write_channel_altera(chan_GA2Conf_active, active);
		write_channel_altera(chan_GA2Conf_mode,   mode);
		write_channel_altera(chan_GA2Conf_cnt,    pop_cnt);
		// --------------------------------------------------------------
		
		#if defined (DEBUG_LEO)
		printf("pop_cnt (INC): %u\n", pop_cnt);
		#endif
	} // End of LOOP_GEN_GENERATIONAL_1
	
	eval_cnt = pop_cnt++;
	
	#if defined (DEBUG_LEO)
	printf("eval_cnt (INC): %u\n", eval_cnt);
	#endif


















	// --------------------------------------------------------------
	// Send DUMMY genotypes to channel to signal INI stop
	// --------------------------------------------------------------
	for (uint pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
		write_channel_altera(chan_GA2Conf_genotype, 0);
	}
	active = 1;
	mode   = 0;
	write_channel_altera(chan_GA2Conf_active, active);
	write_channel_altera(chan_GA2Conf_mode,   mode);
	write_channel_altera(chan_GA2Conf_cnt,    0);

	ack = read_channel_altera(chan_Store2GA_ack);
	//printf("INI ack: %u\n", ack);
	// --------------------------------------------------------------
		











	// Find_best
	uint best_entity_id;
	__local float loc_energies[MAX_POPSIZE];

	// Binary tournament
	uint parent1, parent2;
	__local float local_entity_1     [ACTUAL_GENOTYPE_LENGTH];
	__local float local_entity_2     [ACTUAL_GENOTYPE_LENGTH];	
	__local float offspring_genotype [ACTUAL_GENOTYPE_LENGTH];

	// local search
	uint entity_for_ls;
	uint LS_eval;
	uint evals_for_ls_in_this_cycle;
	uint num_of_evals_for_ls  = 0;
	__local float local_entity_ls	 [ACTUAL_GENOTYPE_LENGTH];
	__local float entity_possible_new_genotype [ACTUAL_GENOTYPE_LENGTH];
	float local_entity_energy;



	// **********************************************
	// ADD VENDOR SPECIFIC PRAGMA
	// **********************************************
	float avg_energy;
	
	LOOP_WHILE_GEN_GENERATIONAL_1:
	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	while ((eval_cnt < DockConst->num_of_energy_evals) && (generation_cnt < DockConst->num_of_generations))
	{
		//Creating a new population

		//Identifying best entity
		best_entity_id = find_best(GlobEnergyCurrent, loc_energies, DockConst->pop_size); 

		#if defined (DEBUG_KERNEL1)
		avg_energy = 0.0f;
		printf("\n\n\nFinal state of the %u. generation:\n", generation_cnt);
		printf("----------------------------\n\n");
		for (i=0; i<DockConst->pop_size; i++)
		{
			avg_energy += GlobEnergyCurrent [i];
			printf("Entity %3u: ", i);
			for (j=0; j<ACTUAL_GENOTYPE_LENGTH; j++) {
				printf("%8.3f ", GlobPopulationCurrent [i*ACTUAL_GENOTYPE_LENGTH +j]);
			}
			printf("   energy sum: %10.3f\n", GlobEnergyCurrent [i*40 + 38]);
		}
		printf("\nAverage energy: %f\nBest energy sum: %f)\n\n", avg_energy/DockConst->pop_size, 
		      GlobEnergyCurrent[best_entity_id]);
		#endif

		//elitism - copying the best entity to new population
		for (uint i=0; i<ACTUAL_GENOTYPE_LENGTH; i++)
		{
			GlobPopulationNext[i] = GlobPopulationCurrent[best_entity_id*ACTUAL_GENOTYPE_LENGTH+i];
		}
		GlobEnergyNext[best_entity_id] = GlobEnergyCurrent[best_entity_id];


		//new population consists of one member currently
		new_pop_cnt = 1;







		// -----------------------------------------------------------------------
		// **********************************************
		// ADD VENDOR SPECIFIC PRAGMA
		// **********************************************
		LOOP_GEN_GENERATIONAL_2:
		for (new_pop_cnt = 1; new_pop_cnt < DockConst->pop_size; new_pop_cnt++)
		{
			//printf("BEFORE BINARY TOURNAMENT, %u\n",new_pop_cnt);
			//selecting two individuals randomly
			binary_tournament_selection(GlobEnergyCurrent,
						    GlobPRNG,
						    &parent1,
						    &parent2,
						    DockConst->pop_size,
						    DockConst->tournament_rate);
			//printf("AFTER BINARY TOURNAMENT, %u\n",new_pop_cnt);


			//mating parents	
			async_work_group_copy(local_entity_1, GlobPopulationCurrent+parent1*ACTUAL_GENOTYPE_LENGTH, ACTUAL_GENOTYPE_LENGTH, 0);
			async_work_group_copy(local_entity_2, GlobPopulationCurrent+parent2*ACTUAL_GENOTYPE_LENGTH, ACTUAL_GENOTYPE_LENGTH, 0);


			//printf("BEFORE GEN NEW GENOTYPE, %u\n",new_pop_cnt);
			// first two args are population [parent1], population [parent2]
			gen_new_genotype(GlobPRNG,
                                         local_entity_1, 
					 local_entity_2,
					 DockConst->mutation_rate,
					 DockConst->abs_max_dmov,
					 DockConst->abs_max_dang,
					 DockConst->crossover_rate,
					 ACTUAL_GENOTYPE_LENGTH,
					 offspring_genotype);

			//printf("AFTER GEN NEW GENOTYPE, %u\n",new_pop_cnt);


			//printf("BEFORE GA CHANNEL, %u\n",new_pop_cnt);
			// --------------------------------------------------------------
			// Send genotypes to channel
			// --------------------------------------------------------------

			for (uint pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
				write_channel_altera(chan_GA2Conf_genotype, offspring_genotype[pipe_cnt]);
			}

			active = 1;
			mode = 2;
			write_channel_altera(chan_GA2Conf_active, active);
			write_channel_altera(chan_GA2Conf_mode,   mode);
			write_channel_altera(chan_GA2Conf_cnt,    new_pop_cnt);

			// --------------------------------------------------------------
			//printf("AFTER GA CHANNEL, %u\n",new_pop_cnt);

			//copying offspring to population
			// **********************************************
			// ADD VENDOR SPECIFIC PRAGMA
			// **********************************************
			
			//LOOP_GEN_GENERATIONAL_3:
			//for (uint gene_cnt = 0; gene_cnt < ACTUAL_GENOTYPE_LENGTH; gene_cnt++) {
			//	GlobPopulationNext[new_pop_cnt*ACTUAL_GENOTYPE_LENGTH + gene_cnt] = offspring_genotype[gene_cnt];
			//}
			
			async_work_group_copy(GlobPopulationNext+new_pop_cnt*ACTUAL_GENOTYPE_LENGTH, offspring_genotype, ACTUAL_GENOTYPE_LENGTH, 0);
	
			#if defined (DEBUG_LEO)
			printf("eval_cnt (INC): %u, new_pop_cnt (INC): %u\n", eval_cnt, new_pop_cnt);
			#endif
			


		
		} // End of LOOP_GEN_GENERATIONAL_2
		// -----------------------------------------------------------------------





		


		#if defined (DEBUG_LEO)
		printf("End of loop of new_pop_cnt, new_pop_cnt = %u\n", new_pop_cnt);
		#endif







		// --------------------------------------------------------------
		// Send DUMMY genotypes to channel to signal GG stop
		// --------------------------------------------------------------
		for (uint pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
			write_channel_altera(chan_GA2Conf_genotype, 0);
		}
		active = 1;
		mode   = 0;
		write_channel_altera(chan_GA2Conf_active, active);
		write_channel_altera(chan_GA2Conf_mode,   mode);
		write_channel_altera(chan_GA2Conf_cnt,    0);

		ack = read_channel_altera(chan_Store2GA_ack);
		//printf("GG ack: %u\n", ack);
		// --------------------------------------------------------------
		




		// Updating old population with new one
		// **********************************************
		// ADD VENDOR SPECIFIC PRAGMA
		// **********************************************
		for (uint i=0;i<DockConst->pop_size*ACTUAL_GENOTYPE_LENGTH; i++)
		{
			GlobPopulationCurrent[i] = GlobPopulationNext[i];
		}

		// Updating old energy with new one
		// **********************************************
		// ADD VENDOR SPECIFIC PRAGMA
		// **********************************************
		for (uint i=0;i<DockConst->pop_size; i++)
		{
			GlobEnergyCurrent[i] = GlobEnergyNext[i];
		}

























		evals_for_ls_in_this_cycle = 0;

		#if defined (DEBUG_LEO)
		printf("num_of_entity_for_ls: %u \n", DockConst->num_of_lsentities);
		#endif









		// -----------------------------------------------------------------------
		// **********************************************
		// ADD VENDOR SPECIFIC PRAGMA
		// **********************************************
		LOOP_GEN_GENERATIONAL_4:
		//subjecting num_of_entity_for_ls pieces of offsprings to LS
		for (uint ls_ent_cnt=0; ls_ent_cnt<DockConst->num_of_lsentities; ls_ent_cnt++)	
		{
			//choosing an entity randomly,
			//and without checking if it has already been subjected to LS in this cycle
			entity_for_ls = myrand_uint(GlobPRNG, DockConst->pop_size);

			#if defined (DEBUG_LEO)
			printf("entity_for_ls: %u\n", entity_for_ls);
			#endif


			#if defined (DEBUG_KERNEL1)
			printf("Entity %u before local search: ", entity_for_ls);
			for (j=0; j<DockConst->rotbondlist_length+6; j++) {
				printf("%f ", GlobPopulationCurrent [entity_for_ls*ACTUAL_GENOTYPE_LENGTH + j]);
			}
			printf("   energies: %f \n", GlobEnergyCurrent [entity_for_ls]);
			#endif

			//printf("BEFORE LS, %u\n",ls_ent_cnt);
			//performing local search
			//async_work_group_copy(local_entity_ls, GlobPopulationCurrent+entity_for_ls*ACTUAL_GENOTYPE_LENGTH, ACTUAL_GENOTYPE_LENGTH, 0);
			//local_entity_energy = GlobEnergyCurrent[entity_for_ls];

			perform_LS(GlobPopulationCurrent,
				   GlobEnergyCurrent,
				   GlobPRNG,
				   KerConst,
				   DockConst,
				   entity_for_ls,
				   local_entity_ls,
				   entity_possible_new_genotype,
				   &LS_eval);

			//async_work_group_copy(GlobPopulationCurrent+entity_for_ls*ACTUAL_GENOTYPE_LENGTH, local_entity_ls, ACTUAL_GENOTYPE_LENGTH, 0);

			//printf("AFTER LS, %u\n",ls_ent_cnt);


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

		//#if defined (DEBUG_LEO)
		printf("eval_cnt: %u, generation_cnt: %u\n", eval_cnt, generation_cnt);
		//#endif




	} // End of while ((eval_cnt < mypars->num_of_energy_evals) && (generation_cnt < mypars->num_of_generations))


	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//*/


















	printf("	%-20s: %s\n", "Krnl_GA", "has finished execution!");

	
	active = 0;
	mode = 4;

	// --------------------------------------------------------------
	// Send DUMMY genotypes to channel
	// --------------------------------------------------------------
	for (uint pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
		write_channel_altera(chan_GA2Conf_genotype, 0);
	}
	write_channel_altera(chan_GA2Conf_active, active);
	write_channel_altera(chan_GA2Conf_mode,   mode);
	write_channel_altera(chan_GA2Conf_cnt,    0);
	// --------------------------------------------------------------


	GlobEvals_performed[0]       =  eval_cnt;
	GlobGenerations_performed[0] =  generation_cnt;

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
#include "Krnl_Store.cl"

