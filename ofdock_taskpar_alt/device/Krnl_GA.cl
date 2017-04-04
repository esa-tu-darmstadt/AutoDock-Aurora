// Enable the channels extension
#pragma OPENCL EXTENSION cl_altera_channels : enable

//IC: initial calculation of energy of populations
//GG: genetic generation 
//LS: local search




// Define kernel file-scope channel variable
// Buffered channels 
// MAX_NUM_OF_ATOMS=90
// ACTUAL_GENOTYPE_LENGTH (MAX_NUM_OF_ROTBONDS+6) =38

channel char  chan_GA2IC_active;
channel uint  chan_IC2GA_eval_cnt;

channel char  chan_GA2GG_active;
channel uint  chan_GG2GA_eval_cnt;

channel char  chan_GA2LS_active;
channel uint  chan_LS2GA_eval_cnt;

// active 1: receiving Kernel is active, 0 receiving Kernel is disabled
// mode 1 for I: init calculation energy, 2 for G: genetic generation, 3 for L: local search
// cnt: population count

channel char  	chan_IC2Conf_active;	
channel char  	chan_IC2Conf_mode;	
channel ushort	chan_IC2Conf_cnt;		
channel float  	chan_IC2Conf_genotype __attribute__((depth(38)));

channel char  	chan_GG2Conf_active;
channel char  	chan_GG2Conf_mode;
channel ushort 	chan_GG2Conf_cnt;
channel float  	chan_GG2Conf_genotype __attribute__((depth(38)));

channel char  	chan_LS2Conf_active;
channel char  	chan_LS2Conf_mode;
channel ushort 	chan_LS2Conf_cnt;
channel float  	chan_LS2Conf_genotype __attribute__((depth(38)));

channel float 	chan_Conf2Intere_x __attribute__((depth(90)));
channel float 	chan_Conf2Intere_y __attribute__((depth(90)));
channel float 	chan_Conf2Intere_z __attribute__((depth(90)));
channel char  	chan_Conf2Intere_active;	
channel char  	chan_Conf2Intere_mode;	
channel ushort  chan_Conf2Intere_cnt;	

channel float 	chan_Conf2Intrae_x __attribute__((depth(90)));
channel float 	chan_Conf2Intrae_y __attribute__((depth(90)));
channel float 	chan_Conf2Intrae_z __attribute__((depth(90)));
channel char  	chan_Conf2Intrae_active;	
channel char  	chan_Conf2Intrae_mode;	
channel ushort  chan_Conf2Intrae_cnt;

channel float 	chan_Intere2Store_intere;
channel char 	chan_Intere2Store_active;	
channel char  	chan_Intere2Store_mode;	
channel ushort  chan_Intere2Store_cnt;

channel float 	chan_Intrae2Store_intrae;
channel char  	chan_Intrae2Store_active;	
channel char  	chan_Intrae2Store_mode;	
channel ushort  chan_Intrae2Store_cnt;

channel char  	chan_Store2IC_ack;
channel char  	chan_Store2GG_ack;
channel float 	chan_Store2LS_LSenergy;

#include "../defines.h"

// Next structures were copied from calcenergy.h

typedef struct
{
	unsigned char  	num_of_atoms;
	unsigned char   num_of_atypes;
	unsigned int    num_of_intraE_contributors;
	unsigned char   gridsize_x;
	unsigned char   gridsize_y;
	unsigned char   gridsize_z;
	unsigned char	g1;
	unsigned int	g2;
	unsigned int 	g3;
	float  		grid_spacing;

	unsigned int    rotbondlist_length;
	float  		coeff_elec;
	float  		coeff_desolv;

	//float* 		conformations_current;
	//float* 		energies_current;
	//float* 		conformations_next;
	//float* 		energies_next;
	//int*   		evals_of_new_entities;
	//unsigned int* 	prng_states;

	// L30nardoSV added
	unsigned int    num_of_energy_evals;
	unsigned int    num_of_generations;

	unsigned int 	pop_size;
	unsigned int	num_of_genes;
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
/*
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
*/
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
} kernelconstant_static;


typedef struct
{
       float ref_coords_x_const[MAX_NUM_OF_ATOMS];
       float ref_coords_y_const[MAX_NUM_OF_ATOMS];
       float ref_coords_z_const[MAX_NUM_OF_ATOMS];
       float rotbonds_moving_vectors_const[3*MAX_NUM_OF_ROTBONDS];
       float rotbonds_unit_vectors_const  [3*MAX_NUM_OF_ROTBONDS];
       float ref_orientation_quats_const  [4];
} kernelconstant_dynamic;

#include "auxiliary_genetic.cl"
//#include "auxiliary_performls.cl"


// --------------------------------------------------------------------------
// The function performs a generational genetic algorithm based search 
// on the search space.
// The first parameter is the population which must be filled with initial values 
// before calling this function. 
// The other parameters are variables which describe the grids, 
// the docking parameters and the ligand to be docked. 
// Originally from: searchoptimum.c
// --------------------------------------------------------------------------
__kernel __attribute__ ((max_global_work_dim(0)))
void Krnl_GA(__global       float*           restrict GlobPopulationCurrent,
	     __global       float*           restrict GlobEnergyCurrent,
	     __global const float*           restrict GlobPopulationNext,
	     __global const float*           restrict GlobEnergyNext,
	     __constant     Dockparameters*  restrict DockConst,	
	     __global       unsigned int*    restrict GlobEvalsGenerations_performed
	     )
{
	//Print algorithm parameters

	#if defined (DEBUG_KRNL_GA)
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

	uint eval_cnt = 0;
	uint generation_cnt = 1;

	char active = 1;

	// IC :  Init Calculation
	write_channel_altera(chan_GA2IC_active, active);
	mem_fence(CLK_CHANNEL_MEM_FENCE);
	eval_cnt = read_channel_altera(chan_IC2GA_eval_cnt);
	mem_fence(CLK_CHANNEL_MEM_FENCE);

	bool GG_act  = false;
	bool GG_done = false;
	bool LS_act  = false;
	bool LS_done = false; 

	uint array_evals_and_generations_performed [2];

	while (active) {
		if ((eval_cnt < DockConst->num_of_energy_evals) && (generation_cnt < DockConst->num_of_generations)) {

			if (LS_done ==  true) {		
				generation_cnt++;

				#if defined (DEBUG_KRNL_GA)
				printf("eval_cnt: %u, generation_cnt: %u\n", eval_cnt, generation_cnt);
				#endif
			}

			//active = 1;
			
			GG_act = true;
			LS_act = false;
			
			if (GG_done == true) {
/*
				// Updating old population with new one
				//for (uint i=0;i<DockConst->pop_size*ACTUAL_GENOTYPE_LENGTH; i++) {
				for (ushort i=0;i<DockConst->pop_size*ACTUAL_GENOTYPE_LENGTH; i++) { 	
					GlobPopulationCurrent[i] = GlobPopulationNext[i];}

				// Updating old energy with new one
				//for (uint i=0;i<DockConst->pop_size; i++) {
				for (ushort i=0;i<DockConst->pop_size; i++) { 			
					GlobEnergyCurrent[i] = GlobEnergyNext[i];}
*/

				for (ushort i=0;i<DockConst->pop_size*ACTUAL_GENOTYPE_LENGTH; i++) { 	
					GlobPopulationCurrent[i] = GlobPopulationNext[i];
					if (i<DockConst->pop_size) {
						GlobEnergyCurrent[i] = GlobEnergyNext[i];
					}
				}

				// LS
				LS_act = true;
			}
		}
		else {
			active = 0; 
		
			#if defined (DEBUG_ACTIVE_KERNEL)
			printf("	%-20s: %s\n", "Krnl_GA", "must be disabled");
			#endif
			
			// Disable Krnl_IC: Automatically
			GG_act  = true;
			GG_done = false;
			LS_act  = true;
			LS_done = false;
			
			array_evals_and_generations_performed [0] = eval_cnt;
			array_evals_and_generations_performed [1] = generation_cnt;
	
			//#pragma unroll 1
			for (uchar i=0; i<2; i++) {
				GlobEvalsGenerations_performed[i] = array_evals_and_generations_performed [i];
			}
		}

		if ((GG_act == true) && (GG_done == false)) {
			write_channel_altera(chan_GA2GG_active, active);
			mem_fence(CLK_CHANNEL_MEM_FENCE);
			eval_cnt += read_channel_altera(chan_GG2GA_eval_cnt);
			mem_fence(CLK_CHANNEL_MEM_FENCE);
			GG_done = true;
			LS_done = false;
		}

		if ((LS_act == true) && (LS_done == false)) {
			write_channel_altera(chan_GA2LS_active, active);
			mem_fence(CLK_CHANNEL_MEM_FENCE);
			eval_cnt += read_channel_altera(chan_LS2GA_eval_cnt);
			mem_fence(CLK_CHANNEL_MEM_FENCE);
			LS_done = true;
			GG_done = false;
		}
		
	} // End while eval_cnt & generation_cnt

	#if defined (DEBUG_ACTIVE_KERNEL)
	printf("	%-20s: %s\n", "Krnl_GA", "disabled");
	#endif
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

#include "Krnl_IC.cl"
#include "Krnl_GG.cl"
#include "Krnl_LS.cl"
#include "Krnl_Conform.cl"
#include "Krnl_InterE.cl"
#include "Krnl_IntraE.cl"
#include "Krnl_Store.cl"

