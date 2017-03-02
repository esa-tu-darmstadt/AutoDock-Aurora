/*
 * (C) 2013. Evopro Innovation Kft.
 *
 * performdocking.c
 *
 *  Created on: 2009.05.25.
 *      Author: pechan.imre
 */

//// ------------------------ 
//// Correct time measurement 
//// Moved from main.cpp to performdocking.cpp 
//// to skip measuring build time 
//#include <sys/time.h> 
//// ------------------------
#include "performdocking.h"

//// -------------------------------- 
//// Altera OpenCL Helper Variables 
//// -------------------------------- 
#include <assert.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <cstring> 
#include "CL/opencl.h" 
#include "AOCLUtils/aocl_utils.h"

using namespace aocl_utils; 

#define STRING_BUFFER_LEN 1024

// OpenCL runtime configuration
static cl_platform_id   platform      = NULL;
static cl_device_id     device        = NULL;
static cl_context       context       = NULL;

// Kernel name, as defined in the CL file
#ifdef ENABLE_KERNEL1
static cl_command_queue command_queue1 = NULL;
static cl_kernel kernel1  = NULL;
static const char *name_k1 = "Krnl_GA";
#endif

#ifdef ENABLE_KERNEL2
static cl_command_queue command_queue2 = NULL;
static cl_kernel kernel2  = NULL;
static const char *name_k2 = "Krnl_Conform";
#endif

#ifdef ENABLE_KERNEL3
static cl_command_queue command_queue3 = NULL;
static cl_kernel kernel3  = NULL;
static const char *name_k3 = "Krnl_InterE";
#endif

#ifdef ENABLE_KERNEL4
static cl_command_queue command_queue4 = NULL;
static cl_kernel kernel4  = NULL;
static const char *name_k4 = "Krnl_IntraE";
#endif

static cl_program program = NULL;

// Function prototypes
bool init();
void cleanup();
static void device_info_ulong  ( cl_device_id device, cl_device_info param, const char* name);
static void device_info_uint   ( cl_device_id device, cl_device_info param, const char* name);
static void device_info_bool   ( cl_device_id device, cl_device_info param, const char* name);
static void device_info_string ( cl_device_id device, cl_device_info param, const char* name);
static void display_device_info( cl_device_id device );

//// --------------------------------
//// Host constant struct
//// --------------------------------
Dockingconstant DockConst;
Ligandconstant  LigConst;
Gridconstant    GridConst;

//// --------------------------------
//// Host memory buffers
//// --------------------------------


//// --------------------------------
//// Device memory buffers
//// --------------------------------
cl_mem mem_DockConst;			
cl_mem mem_LigConst;
cl_mem mem_GridConst;

// Krnl_GA
cl_mem mem_GlobPopulation;
cl_mem mem_GlobLigandAtom_idxyzq;
cl_mem mem_GlobFgrids;
cl_mem mem_GlobCounter;
cl_mem mem_GlobPRNG;

// Krnl_Conform
cl_mem mem_GlobLigand_rotbonds;
cl_mem mem_GlobLigand_atom_rotbonds;
cl_mem mem_GlobLigand_rotbonds_moving_vectors;
cl_mem mem_GlobLigand_rotbonds_uint_vectors;

// Krnl_InterE
//cl_mem mem_GlobFgrids;

// Krnl_IntraE
cl_mem mem_GlobLigand_intraE_contributors;
cl_mem mem_GlobLigand_VWpars_A;
cl_mem mem_GlobLigand_VWpars_B;
cl_mem mem_GlobLigand_VWpars_C;
cl_mem mem_GlobLigand_VWpars_D;
cl_mem mem_GlobLigand_VWpars_volume;
cl_mem mem_GlobLigand_VWpars_solpar;
cl_mem mem_Glob_r_6_table;
cl_mem mem_Glob_r_10_table;
cl_mem mem_Glob_r_12_table;
cl_mem mem_Glob_r_epsr_table;
cl_mem mem_Glob_desolv_table;
cl_mem mem_Glob_q1q2;
cl_mem mem_Glob_qasp_mul_absq;

//// --------------------------------
//// Docking
//// --------------------------------
/*
int docking_with_cpu(const Gridinfo* mygrid, const double* floatgrids, Dockpars* mypars, const Liganddata* myligand_init,
					  const int* argc, char** argv, double* cpu_exectime_sum, Ligandresult* result_ligands)
*/

int docking_with_cpu(const Gridinfo* mygrid, 		     
		     const float* floatgrids, 		     
		     Dockpars* mypars, 		     
		     const Liganddata* myligand_init, 		     
		     const int* argc, 		    
		     char** argv, 		     
		     double* cpu_exectime_sum,	// this can be double 		     
		     Ligandresult* result_ligands)


//The function performs the docking algorithm with the CPU.
{
	////------------------------------------------------------------ 	
	//// OpenCL Host Setup         
	////------------------------------------------------------------ 	
	if(!init()) {     		
		return -1;   	
	} 	
	printf("Init complete!\n"); fflush(stdout);

	// Check OpenCL call exhaustively
	cl_int status;
	////------------------------------------------------------------ 




	double start_algorithm, stop_algorithm;
	int run_cnt;
	int i, j;

	Liganddata myligand;

	float movvec_to_origo [3];
	float init_population [CPU_MAX_POP_SIZE][40];
	float final_population [CPU_MAX_POP_SIZE][40];
	unsigned int termination_counter [2] = {0, 0};


	*cpu_exectime_sum = 0;

	int debug = 0;

	// Initial load, required to set scalar values of KERNEL1 	
	myligand = *myligand_init;

	// PRNGs values need to be generated 	
	unsigned int prng[40]; 	
	srand((unsigned int) time(NULL)); 	
	for (i=0; i<40;i++) { 
#if defined (REPRO) 		
		prng[i] = 1u; 
#else 		
		prng[i] = rand(); 		
		//printf("prng[%u]: %u\n", i, prng[i]); 
#endif 	
	}


	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// to simplify intraE computation in kernel
	// extracted from ofdock_optimized_3/calcenergy.cpp
	int type_id1, type_id2;
	char* charpoi;
	char  intraE_contributors[3*MAX_INTRAE_CONTRIBUTORS];

	//intramolecular energy contributors
	myligand.num_of_intraE_contributors = 0;
	for (i=0; i<myligand.num_of_atoms-1; i++)
		for (j=i+1; j<myligand.num_of_atoms; j++)
		{
			if (myligand.intraE_contributors[i][j])
				myligand.num_of_intraE_contributors++;
		}

	if (myligand.num_of_intraE_contributors > MAX_INTRAE_CONTRIBUTORS)
	{
		printf("Error: number of intramolecular energy contributor is too high!\n");
		return 1;
	}

	charpoi = intraE_contributors;
	for (i=0; i<myligand.num_of_atoms-1; i++)
		for (j=i+1; j<myligand.num_of_atoms; j++)
		{
			if (myligand.intraE_contributors[i][j] == 1)
			{
				*charpoi = (char) i;
				charpoi++;
				*charpoi = (char) j;
				charpoi++;

				type_id1 = (int) myligand.atom_idxyzq [i][0];
				type_id2 = (int) myligand.atom_idxyzq [j][0];
				if (is_H_bond(myligand.atom_types[type_id1], myligand.atom_types[type_id2]) != 0)
					*charpoi = (char) 1;
				else
					*charpoi = (char) 0;
				charpoi++;
			}
		}

	//------------------------------------------------------------
	float r_6_table    [2048];
	float r_10_table   [2048];
	float r_12_table   [2048];
	float r_epsr_table [2048];
	float desolv_table [2048];

	float dist = 0.0f;
	float calc_ddd_Mehler_Solmajer;
	
	for (i=0; i<2048; i++) {
		dist += 0.01f;
		r_6_table  [i] = 1/pow(dist,6);
		r_10_table [i] = 1/pow(dist,10);
		r_12_table [i] = 1/pow(dist,12);
		calc_ddd_Mehler_Solmajer = A + (B / (1.0f + RK*exp(LAMBDA_B * dist)));
		r_epsr_table [i] = mypars->coeffs.scaled_AD4_coeff_elec / (dist*calc_ddd_Mehler_Solmajer);
		desolv_table [i] = mypars->coeffs.AD4_coeff_desolv*exp( (-1*pow(dist,2)) / (2*pow(SIGMA,2)) );
	}

	//------------------------------------------------------------
	float q1q2 [MAX_NUM_OF_ATOMS*MAX_NUM_OF_ATOMS];

	for (i=0; i<MAX_NUM_OF_ATOMS; i++) {
		for (j=0; j<MAX_NUM_OF_ATOMS; j++) {
			q1q2 [i*MAX_NUM_OF_ATOMS+j] = myligand.atom_idxyzq[i][4] * myligand.atom_idxyzq[j][4];
		}
	}

	//------------------------------------------------------------
	float qasp_mul_absq [MAX_NUM_OF_ATOMS];

	for (i=0; i<MAX_NUM_OF_ATOMS; i++) {
		qasp_mul_absq [i] = mypars->qasp*fabs(myligand.atom_idxyzq[i][4]);
	}

	//------------------------------------------------------------
/*
	// to extract charges and pass them before main loop 
	float atom_charges[MAX_NUM_OF_ATOMS];

	for (i=0; i < myligand.num_of_atoms; i++)
	{
		atom_charges[i] = (float) myligand.atom_idxyzq[i][4];
	}
*/

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	


	DockConst.num_of_energy_evals = mypars->num_of_energy_evals;
	DockConst.num_of_generations  = mypars->num_of_generations;
	DockConst.dmov_mask 	      = mypars->dmov_mask / pow(2.0f, 10.0f);
	DockConst.dang_mask 	      = mypars->dang_mask / pow(2.0f, 8.0f) * 180.0f/ 512.0f;
	DockConst.mutation_rate       = mypars->mutation_rate / 255.0f;
	DockConst.crossover_rate      = mypars->crossover_rate / 255.0f;
	DockConst.lsearch_rate 	      = mypars->lsearch_rate;
	DockConst.tournament_rate     = mypars->tournament_rate / 255.0f;
	DockConst.rho_lower_bound     = mypars->rho_lower_bound / pow(2.0f, 10.0f);
	DockConst.base_dmov_mul_sqrt3 = mypars->base_dmov_mul_sqrt3 / pow(2.0f, 10.0f);
	DockConst.base_dang_mul_sqrt3 = mypars->base_dang_mul_sqrt3 * 180.0f/ 512.0f/ pow(2.0f, 8.0f);
        DockConst.cons_limit 	      = mypars->cons_limit + 1;
        DockConst.max_num_of_iters    = mypars->max_num_of_iters;
	DockConst.pop_size 	      = mypars->pop_size + 1;
	DockConst.num_of_entity_for_ls= (unsigned int) floor((DockConst.pop_size + 1) * DockConst.lsearch_rate + 0.5f);

	LigConst.num_of_atoms 	            = myligand.num_of_atoms;
	LigConst.num_of_rotbonds 	    = myligand.num_of_rotbonds;
	LigConst.num_of_intraE_contributors = myligand.num_of_intraE_contributors;

	GridConst.size_x 		    = mygrid->size_xyz[0];
	GridConst.size_y 		    = mygrid->size_xyz[1];
	GridConst.size_z 		    = mygrid->size_xyz[2];
	GridConst.g1 		            = mygrid->size_xyz[0];
	GridConst.g2 		            = mygrid->size_xyz[0] * mygrid->size_xyz[1];
	GridConst.g3 		    	    = mygrid->size_xyz[0] * mygrid->size_xyz[1] * mygrid->size_xyz[2];
	GridConst.spacing 	            = mygrid->spacing;
	GridConst.num_of_atypes 	    = mygrid->num_of_atypes;


	printf("GridConst.size_x = %u\n", GridConst.size_x);	
	printf("GridConst.size_y = %u\n", GridConst.size_y);
	printf("GridConst.size_z = %u\n", GridConst.size_z);
	printf("GridConst.g1 = %u\n", GridConst.g1);
	printf("GridConst.g2 = %u\n", GridConst.g2);
	printf("GridConst.g3 = %u\n", GridConst.g3);
	printf("GridConst.spacing  = %f\n", GridConst.spacing);
	printf("GridConst.num_of_atypes = %u\n", GridConst.num_of_atypes);

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//------------------------------------------------------------
	// Krnl_GA
	//------------------------------------------------------------

	size_t size_population = CPU_MAX_POP_SIZE*40*sizeof(float);
	size_t size_ligandatom_idxyzq = MAX_NUM_OF_ATOMS*5*sizeof(float);
	size_t size_floatgrids = (mygrid->num_of_atypes+2)*(mygrid->size_xyz[0])*(mygrid->size_xyz[1])*(mygrid->size_xyz[2])*sizeof(float);
	size_t size_counter    = 2*sizeof(unsigned int);
	size_t size_prng       = 40*sizeof(int);

	status = mallocBufferObject(context,CL_MEM_READ_WRITE, size_population,        &mem_GlobPopulation);
	checkError(status, "Failed to allocate: mem_GlobPopulation\n");
	status = mallocBufferObject(context,CL_MEM_READ_ONLY,  size_ligandatom_idxyzq, &mem_GlobLigandAtom_idxyzq);
	checkError(status, "Failed to allocate: mem_GlobLigandAtom_idxyzq\n");
	status = mallocBufferObject(context,CL_MEM_READ_ONLY,  size_floatgrids,        &mem_GlobFgrids);
	checkError(status, "Failed to allocate: mem_GlobFgrids\n");
	status = mallocBufferObject(context,CL_MEM_READ_ONLY,  size_counter, 	       &mem_GlobCounter);
	checkError(status, "Failed to allocate: mem_GlobCounter\n");
	status = mallocBufferObject(context,CL_MEM_READ_ONLY,  size_prng,	       &mem_GlobPRNG);
	checkError(status, "Failed to allocate: mem_GlobPRNG\n");
	status = mallocBufferObject(context,CL_MEM_READ_ONLY,  sizeof(DockConst),      &mem_DockConst);
	checkError(status, "Failed to allocate: mem_DockConst\n");
	status = mallocBufferObject(context,CL_MEM_READ_ONLY,  sizeof(LigConst),       &mem_LigConst);
	checkError(status, "Failed to allocate: mem_LigConst\n");
	status = mallocBufferObject(context,CL_MEM_READ_ONLY,  sizeof(GridConst),      &mem_GridConst);
	checkError(status, "Failed to allocate: mem_GridConst\n");
#ifdef ENABLE_KERNEL1
	status = memcopyBufferObjectToDevice(command_queue1,mem_GlobPopulation,        final_population,     size_population);
	checkError(status, "Failed to copy data to device: mem_GlobPopulation\n");
	status = memcopyBufferObjectToDevice(command_queue1,mem_GlobLigandAtom_idxyzq, myligand.atom_idxyzq, size_ligandatom_idxyzq);
	checkError(status, "Failed to copy data to device: mem_GlobLigandAtom_idxyzq\n");
	status = memcopyBufferObjectToDevice(command_queue1,mem_GlobFgrids, 	       (void*) floatgrids,   size_floatgrids);
	checkError(status, "Failed to copy data to device: mem_GlobFgrids\n");
	status = memcopyBufferObjectToDevice(command_queue1,mem_GlobCounter, 	       termination_counter,  size_counter);
	checkError(status, "Failed to copy data to device: mem_GlobCounter\n");
	status = memcopyBufferObjectToDevice(command_queue1,mem_GlobPRNG,              prng,                 size_prng);
	checkError(status, "Failed to copy data to device: mem_GlobPRNG\n");
	status = memcopyBufferObjectToDevice(command_queue1,mem_DockConst,             &DockConst,           sizeof(DockConst));
	checkError(status, "Failed to copy data to device: mem_DockConst\n");
	status = memcopyBufferObjectToDevice(command_queue1,mem_LigConst,              &LigConst,            sizeof(LigConst));
	checkError(status, "Failed to copy data to device: mem_LigConst\n");
	status = memcopyBufferObjectToDevice(command_queue1,mem_GridConst,             &GridConst,           sizeof(GridConst));
	checkError(status, "Failed to copy data to device: mem_GridConst\n");

	status = setKernelArg(kernel1, 0, sizeof(cl_mem), &mem_GlobPopulation);
	checkError(status, "Failed to set krnl1 arg 0\n");
	status = setKernelArg(kernel1, 1, sizeof(cl_mem), &mem_GlobLigandAtom_idxyzq);
	checkError(status, "Failed to set krnl1 arg 1\n");
	status = setKernelArg(kernel1, 2, sizeof(cl_mem), &mem_GlobFgrids);
	checkError(status, "Failed to set krnl1 arg 2\n");
	status = setKernelArg(kernel1, 3, sizeof(cl_mem), &mem_GlobCounter);
	checkError(status, "Failed to set krnl1 arg 3\n");
	status = setKernelArg(kernel1, 4, sizeof(cl_mem), &mem_GlobPRNG);
	checkError(status, "Failed to set krnl1 arg 4\n");
	status = setKernelArg(kernel1, 5, sizeof(cl_mem), &mem_DockConst);
	checkError(status, "Failed to set krnl1 arg 5\n");
	status = setKernelArg(kernel1, 6, sizeof(cl_mem), &mem_LigConst);
	checkError(status, "Failed to set krnl1 arg 6\n");
	status = setKernelArg(kernel1, 7, sizeof(cl_mem), &mem_GridConst);
	checkError(status, "Failed to set krnl1 arg 7\n");
#endif

	//------------------------------------------------------------
	// Krnl_Conform
	//------------------------------------------------------------	
	size_t size_ligand_rotbonds         = MAX_NUM_OF_ROTBONDS*2*sizeof(int);
	size_t size_ligand_atom_rotbonds    = MAX_NUM_OF_ATOMS*32*sizeof(char);
	size_t size_ligand_rotbonds_vectors = MAX_NUM_OF_ROTBONDS*3*sizeof(float);

	status = mallocBufferObject(context,CL_MEM_READ_WRITE, size_ligand_rotbonds,         &mem_GlobLigand_rotbonds);
	checkError(status, "Failed to allocate: mem_GlobLigand_rotbonds\n");
	status = mallocBufferObject(context,CL_MEM_READ_WRITE, size_ligand_atom_rotbonds,    &mem_GlobLigand_atom_rotbonds);
	checkError(status, "Failed to allocate: mem_GlobLigand_atom_rotbonds\n");
	status = mallocBufferObject(context,CL_MEM_READ_WRITE, size_ligand_rotbonds_vectors, &mem_GlobLigand_rotbonds_moving_vectors);
	checkError(status, "Failed to allocate: mem_GlobLigand_rotbonds_moving_vectors\n");
	status = mallocBufferObject(context,CL_MEM_READ_WRITE, size_ligand_rotbonds_vectors, &mem_GlobLigand_rotbonds_uint_vectors);
	checkError(status, "Failed to allocate: mem_GlobLigand_rotbonds_uint_vectors\n");
#ifdef ENABLE_KERNEL2	
	status = memcopyBufferObjectToDevice(command_queue2,mem_GlobLigand_rotbonds, myligand.rotbonds, size_ligand_rotbonds);
	checkError(status, "Failed to copy data to device: mem_GlobLigand_rotbonds\n");
	status = memcopyBufferObjectToDevice(command_queue2,mem_GlobLigand_atom_rotbonds, myligand.atom_rotbonds, size_ligand_atom_rotbonds);
	checkError(status, "Failed to copy data to device: mem_GlobLigand_atom_rotbonds\n");
	status = memcopyBufferObjectToDevice(command_queue2,mem_GlobLigand_rotbonds_moving_vectors, myligand.rotbonds_moving_vectors, size_ligand_rotbonds_vectors);
	checkError(status, "Failed to copy data to device: mem_GlobLigand_rotbonds_moving_vectors\n");
	status = memcopyBufferObjectToDevice(command_queue2,mem_GlobLigand_rotbonds_uint_vectors, myligand.rotbonds_unit_vectors, size_ligand_rotbonds_vectors);
	checkError(status, "Failed to copy data to device: mem_GlobLigand_rotbonds_uint_vectors\n");

	status = setKernelArg(kernel2, 0, sizeof(cl_mem), &mem_GlobLigand_rotbonds);
	checkError(status, "Failed to set krnl2 arg 0\n");
	status = setKernelArg(kernel2, 1, sizeof(cl_mem), &mem_GlobLigand_atom_rotbonds);
	checkError(status, "Failed to set krnl2 arg 1\n");
	status = setKernelArg(kernel2, 2, sizeof(cl_mem), &mem_GlobLigand_rotbonds_moving_vectors);
	checkError(status, "Failed to set krnl2 arg 2\n");
	status = setKernelArg(kernel2, 3, sizeof(cl_mem), &mem_GlobLigand_rotbonds_uint_vectors);
	checkError(status, "Failed to set krnl2 arg 3\n");
	status = setKernelArg(kernel2, 4, sizeof(cl_mem), &mem_LigConst);
	checkError(status, "Failed to set krnl2 arg 4\n");
#endif

	//------------------------------------------------------------
	// Krnl_InterE
	//------------------------------------------------------------
#ifdef ENABLE_KERNEL3
	setKernelArg(kernel3, 0, sizeof(cl_mem), &mem_GlobFgrids);
	setKernelArg(kernel3, 1, sizeof(cl_mem), &mem_LigConst);
	setKernelArg(kernel3, 2, sizeof(cl_mem), &mem_GridConst);
#endif

	//------------------------------------------------------------
	// Krnl_IntraE
	//------------------------------------------------------------
	size_t size_ligand_intrae_contrib = 3*MAX_INTRAE_CONTRIBUTORS*sizeof(char);
	size_t size_ligand_vwpars         = MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES*sizeof(float);
	size_t size_ligand_vwpars2        = MAX_NUM_OF_ATYPES*sizeof(float);
	size_t size_r_table		  = 2048*sizeof(float);
	size_t size_q1q2		  = MAX_NUM_OF_ATOMS*MAX_NUM_OF_ATOMS*sizeof(float);
	size_t size_qasp_mul_absq         = MAX_NUM_OF_ATOMS*sizeof(float);

	mallocBufferObject(context,CL_MEM_READ_ONLY, size_ligand_intrae_contrib, &mem_GlobLigand_intraE_contributors);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_ligand_vwpars,  &mem_GlobLigand_VWpars_A);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_ligand_vwpars,  &mem_GlobLigand_VWpars_B);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_ligand_vwpars,  &mem_GlobLigand_VWpars_C);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_ligand_vwpars,  &mem_GlobLigand_VWpars_D);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_ligand_vwpars2, &mem_GlobLigand_VWpars_volume);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_ligand_vwpars2, &mem_GlobLigand_VWpars_solpar);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_r_table,        &mem_Glob_r_6_table);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_r_table,        &mem_Glob_r_10_table);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_r_table,        &mem_Glob_r_12_table);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_r_table,        &mem_Glob_r_epsr_table);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_r_table,        &mem_Glob_desolv_table);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_q1q2,           &mem_Glob_q1q2);
	mallocBufferObject(context,CL_MEM_READ_ONLY, size_qasp_mul_absq,  &mem_Glob_qasp_mul_absq);
#ifdef ENABLE_KERNEL4	
	memcopyBufferObjectToDevice(command_queue4,mem_GlobLigand_intraE_contributors, intraE_contributors, size_ligand_intrae_contrib);
	memcopyBufferObjectToDevice(command_queue4,mem_GlobLigand_VWpars_A,      myligand.VWpars_A, size_ligand_vwpars);
	memcopyBufferObjectToDevice(command_queue4,mem_GlobLigand_VWpars_B,      myligand.VWpars_B, size_ligand_vwpars);
	memcopyBufferObjectToDevice(command_queue4,mem_GlobLigand_VWpars_C,      myligand.VWpars_C, size_ligand_vwpars);
	memcopyBufferObjectToDevice(command_queue4,mem_GlobLigand_VWpars_D,      myligand.VWpars_D, size_ligand_vwpars);
	memcopyBufferObjectToDevice(command_queue4,mem_GlobLigand_VWpars_volume, myligand.volume,   size_ligand_vwpars2);
	memcopyBufferObjectToDevice(command_queue4,mem_GlobLigand_VWpars_solpar, myligand.solpar,   size_ligand_vwpars2);
	memcopyBufferObjectToDevice(command_queue4,mem_Glob_r_6_table,           r_6_table,         size_r_table);
	memcopyBufferObjectToDevice(command_queue4,mem_Glob_r_10_table,          r_10_table,        size_r_table);
	memcopyBufferObjectToDevice(command_queue4,mem_Glob_r_12_table,          r_12_table,        size_r_table);
	memcopyBufferObjectToDevice(command_queue4,mem_Glob_r_epsr_table,        r_epsr_table,      size_r_table);
	memcopyBufferObjectToDevice(command_queue4,mem_Glob_desolv_table,        desolv_table,      size_r_table);
	memcopyBufferObjectToDevice(command_queue4,mem_Glob_q1q2, 		 q1q2, 		    size_q1q2);
	memcopyBufferObjectToDevice(command_queue4,mem_Glob_qasp_mul_absq,       qasp_mul_absq,     size_qasp_mul_absq);

	setKernelArg(kernel4, 0, sizeof(cl_mem), &mem_GlobLigand_intraE_contributors);  
	setKernelArg(kernel4, 1, sizeof(cl_mem), &mem_GlobLigand_VWpars_A);
	setKernelArg(kernel4, 2, sizeof(cl_mem), &mem_GlobLigand_VWpars_B);
	setKernelArg(kernel4, 3, sizeof(cl_mem), &mem_GlobLigand_VWpars_C);
	setKernelArg(kernel4, 4, sizeof(cl_mem), &mem_GlobLigand_VWpars_D);
	setKernelArg(kernel4, 5, sizeof(cl_mem), &mem_GlobLigand_VWpars_volume);
	setKernelArg(kernel4, 6, sizeof(cl_mem), &mem_GlobLigand_VWpars_solpar);
	setKernelArg(kernel4, 7, sizeof(cl_mem), &mem_Glob_r_6_table);
	setKernelArg(kernel4, 8, sizeof(cl_mem), &mem_Glob_r_10_table);
	setKernelArg(kernel4, 9, sizeof(cl_mem), &mem_Glob_r_12_table);
	setKernelArg(kernel4,10, sizeof(cl_mem), &mem_Glob_r_epsr_table);
	setKernelArg(kernel4,11, sizeof(cl_mem), &mem_Glob_desolv_table);
	setKernelArg(kernel4,12, sizeof(cl_mem), &mem_Glob_q1q2);
	setKernelArg(kernel4,13, sizeof(cl_mem), &mem_Glob_qasp_mul_absq);
	setKernelArg(kernel4,14, sizeof(cl_mem), &mem_LigConst);
	setKernelArg(kernel4,15, sizeof(cl_mem), &mem_GridConst);
#endif

	////------------------------------------------------------------ 	
	//// Docking 	
	////------------------------------------------------------------
	printf("Starting dockings...\n");
	fflush(stdout);

	for (run_cnt = 0; run_cnt < mypars->num_of_runs; run_cnt++)
	{

		myligand = *myligand_init;

		// Capturing seeds and initial population (getparameters_pipeline.cpp)
			// Writen: mypars: mypars->ref_ori_angles
			//		   mypars->seed
			// Written: init_population
		get_seeds_and_initpop(mypars, init_population, &myligand, mygrid);	

		// Capturing reference orientation of the ligand
			// Writte: myligand: myligand->atom_idxyzq    
		get_ref_orientation(&myligand, mypars);	

		//x, y and z in grid spacing instead of A
		for (i=0; i<CPU_MAX_POP_SIZE; i++)
			for (j=0; j<3; j++)
				init_population [i][j] = init_population [i][j]/mygrid->spacing;	

		//Moving, rotating and scaling the ligand
		get_movvec_to_origo(&myligand, movvec_to_origo);	// Written: movvec_to_origo
		move_ligand(&myligand, movvec_to_origo);		// Written: myligand: myligand->atom_idxyzq
		scale_ligand(&myligand, (double) 1.0/mygrid->spacing);	// Written: myligand: myligand->atom_idxyzq

		//Unit and moving vectors of the ligand's rotatable bonds 
		//must be re-calculated due to the new orientation and scaling
		get_moving_and_unit_vectors(&myligand);			// Written: myligand: myligand->rotbonds_moving_vectors
									//                    myligand->rotbonds_unit_vectors

		//Executing docking algorithm
		memcpy(final_population, init_population, sizeof(double)*CPU_MAX_POP_SIZE*40);


#if defined (DEBUG_HOST)

#endif

		//------------------------------------------------------------
		//Copy Data to Kernel 0 (NEEDED AS SOME DATA WAS JUST GENERATED WITHIN THIS LOOP)
		//------------------------------------------------------------		
#ifdef ENABLE_KERNEL1
		memcopyBufferObjectToDevice(command_queue1,mem_GlobPopulation, final_population, size_population);
	        memcopyBufferObjectToDevice(command_queue1,mem_GlobLigandAtom_idxyzq, myligand.atom_idxyzq, size_ligandatom_idxyzq);
	        memcopyBufferObjectToDevice(command_queue1,mem_GlobLigand_rotbonds_moving_vectors, myligand.rotbonds_moving_vectors, size_ligand_rotbonds_vectors);
	        memcopyBufferObjectToDevice(command_queue1,mem_GlobLigand_rotbonds_uint_vectors,   myligand.rotbonds_unit_vectors,   size_ligand_rotbonds_vectors);
#endif


		//------------------------------------------------------------
		//Enqueueing Kernels
		//------------------------------------------------------------	
		start_algorithm = timer_gets();
		printf("Run %d started...     \n", run_cnt+1);
		fflush(stdout);
#ifdef ENABLE_KERNEL1
	 	runKernelTask(command_queue1,kernel1,NULL,NULL);
#endif

#ifdef ENABLE_KERNEL2
		runKernelTask(command_queue2,kernel2,NULL,NULL);
#endif

#ifdef ENABLE_KERNEL3
		runKernelTask(command_queue3,kernel3,NULL,NULL);
#endif

#ifdef ENABLE_KERNEL4
		runKernelTask(command_queue4,kernel4,NULL,NULL);
#endif

#ifdef ENABLE_KERNEL1
		clFinish(command_queue1);
#endif

#ifdef ENABLE_KERNEL2
		clFinish(command_queue2);
#endif

#ifdef ENABLE_KERNEL3	
		clFinish(command_queue3);
#endif

#ifdef ENABLE_KERNEL4
		clFinish(command_queue4);
#endif

		//------------------------------------------------------------
		//Copy Data from Kernel 2 (NEEDED AS SOME DATA WAS UPDATED WITHIN KERNELS)
		//------------------------------------------------------------
#ifdef ENABLE_KERNEL1
		memcopyBufferObjectFromDevice(command_queue1,final_population,   mem_GlobPopulation,size_population);
		memcopyBufferObjectFromDevice(command_queue1,termination_counter,mem_GlobCounter,   size_counter);
#endif

#if defined (DEBUG_HOST)
		printf("Number of energy evals performed: %u\n", termination_counter[0]);
		printf("Number of generations performed : %u\n", termination_counter[1]);
#endif





		stop_algorithm = timer_gets();

		printf("finished, CPU run time: %.3fs\n", stop_algorithm - start_algorithm);
		fflush(stdout);

		*cpu_exectime_sum += (stop_algorithm - start_algorithm);


		//------------------------------------------------------------
		//Acquiring and processing result
		//------------------------------------------------------------
		arrange_result(final_population, (int) mypars->pop_size + 1);


		make_resfiles(final_population, &myligand, myligand_init, mypars, mygrid, floatgrids,
					  argc, argv, debug, run_cnt, &(result_ligands [run_cnt]));


	}

  	// Free the resources allocated   	
	cleanup();

	return 0;
}

//// --------------------------------
//// Altera OpenCL Helper Functions
//// --------------------------------
bool init() {
  cl_int status;

  if(!setCwdToExeDir()) {
    return false;
  }

  // Get the OpenCL platform.
  platform = findPlatform("Intel(R) FPGA");
  //platform = findPlatform("Altera SDK");
  if(platform == NULL) {
    printf("ERROR: Unable to find Intel(R) FPGA OpenCL platform.\n");
    return false;
  }

  // User-visible output - Platform information
  {
    char char_buffer[STRING_BUFFER_LEN]; 
    printf("Querying platform for info:\n");
    printf("==========================\n");
    clGetPlatformInfo(platform, CL_PLATFORM_NAME, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n", "CL_PLATFORM_NAME", char_buffer);
    clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n", "CL_PLATFORM_VENDOR ", char_buffer);
    clGetPlatformInfo(platform, CL_PLATFORM_VERSION, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n\n", "CL_PLATFORM_VERSION ", char_buffer);
  }

  // Query the available OpenCL devices.
  scoped_array<cl_device_id> devices;
  cl_uint num_devices;

  devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));

  // We'll just use the first device.
  device = devices[0];

  // Display some device information.
  display_device_info(device);

  // Create the context.
  context = clCreateContext(NULL, 1, &device, &oclContextCallback, NULL, &status);
  checkError(status, "Failed to create context");

  // Create the command queue.
  //queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
#ifdef ENABLE_KERNEL1
  command_queue1 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue 1");
#endif

#ifdef ENABLE_KERNEL2
  command_queue2 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue 2");
#endif

#ifdef ENABLE_KERNEL3
  command_queue3 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue 3");
#endif

#ifdef ENABLE_KERNEL4
  command_queue4 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue 4");
#endif

  // Create the program.
  std::string binary_file = getBoardBinaryFile("docking", device);
  printf("Using AOCX: %s\n", binary_file.c_str());
  program = createProgramFromBinary(context, binary_file.c_str(), &device, 1);


  // Build the program that was just created.

  status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
  checkError(status, "Failed to build program");


  // Create the kernel - name passed in here must match kernel name in the
  // original CL file, that was compiled into an AOCX file using the AOC tool
#ifdef ENABLE_KERNEL1
  kernel1 = clCreateKernel(program, name_k1, &status);
  checkError(status, "Failed to create kernel");
#endif

#ifdef ENABLE_KERNEL2
  kernel2 = clCreateKernel(program, name_k2, &status);
  checkError(status, "Failed to create kernel");
#endif

#ifdef ENABLE_KERNEL4
  kernel4 = clCreateKernel(program, name_k4, &status);
  checkError(status, "Failed to create kernel");
#endif

#ifdef ENABLE_KERNEL3
  kernel3 = clCreateKernel(program, name_k3, &status);
  checkError(status, "Failed to create kernel");
#endif

  return true;
}

// Free the resources allocated during initialization
void cleanup() {
#ifdef ENABLE_KERNEL1
  if(kernel1) {clReleaseKernel(kernel1);}
  if(command_queue1) {clReleaseCommandQueue(command_queue1);}
#endif

#ifdef ENABLE_KERNEL2
  if(kernel2) {clReleaseKernel(kernel2);}
  if(command_queue2) {clReleaseCommandQueue(command_queue2);}
#endif

#ifdef ENABLE_KERNEL3
  if(kernel3) {clReleaseKernel(kernel3);}
  if(command_queue3) {clReleaseCommandQueue(command_queue3);}
#endif

#ifdef ENABLE_KERNEL4
  if(kernel4) {clReleaseKernel(kernel4);}
  if(command_queue4) {clReleaseCommandQueue(command_queue4);}
#endif

  if(program) {clReleaseProgram(program);}
  if(context) {clReleaseContext(context);}

  if(mem_DockConst)			{clReleaseMemObject(mem_DockConst);}
  if(mem_LigConst)			{clReleaseMemObject(mem_LigConst);}
  if(mem_GridConst)			{clReleaseMemObject(mem_GridConst);}
  // Krnl_GA
  if(mem_GlobPopulation) 		{clReleaseMemObject(mem_GlobPopulation);}
  if(mem_GlobLigandAtom_idxyzq) 	{clReleaseMemObject(mem_GlobLigandAtom_idxyzq);}
  if(mem_GlobFgrids) 			{clReleaseMemObject(mem_GlobFgrids);}
  if(mem_GlobCounter) 			{clReleaseMemObject(mem_GlobCounter);}
  if(mem_GlobPRNG) 			{clReleaseMemObject(mem_GlobPRNG);}
  // Krnl_Conform
  if(mem_GlobLigand_rotbonds) 			{clReleaseMemObject(mem_GlobLigand_rotbonds);}
  if(mem_GlobLigand_atom_rotbonds) 		{clReleaseMemObject(mem_GlobLigand_atom_rotbonds);}
  if(mem_GlobLigand_rotbonds_moving_vectors) 	{clReleaseMemObject(mem_GlobLigand_rotbonds_moving_vectors);}
  if(mem_GlobLigand_rotbonds_uint_vectors) 	{clReleaseMemObject(mem_GlobLigand_rotbonds_uint_vectors);}
  // Krnl_InterE
  //if(mem_GlobFgrids) 			{clReleaseMemObject(mem_GlobFgrids);}
  // Krnl_IntraE
  if(mem_GlobLigand_intraE_contributors){clReleaseMemObject(mem_GlobLigand_intraE_contributors);}
  if(mem_GlobLigand_VWpars_A) 		{clReleaseMemObject(mem_GlobLigand_VWpars_A);}
  if(mem_GlobLigand_VWpars_B)   	{clReleaseMemObject(mem_GlobLigand_VWpars_B);}
  if(mem_GlobLigand_VWpars_C)        	{clReleaseMemObject(mem_GlobLigand_VWpars_C);}
  if(mem_GlobLigand_VWpars_D)		{clReleaseMemObject(mem_GlobLigand_VWpars_D);}
  if(mem_GlobLigand_VWpars_volume)      {clReleaseMemObject(mem_GlobLigand_VWpars_volume);}
  if(mem_GlobLigand_VWpars_solpar)      {clReleaseMemObject(mem_GlobLigand_VWpars_solpar);}
  if(mem_Glob_r_6_table)      		{clReleaseMemObject(mem_Glob_r_6_table);}
  if(mem_Glob_r_10_table)     		{clReleaseMemObject(mem_Glob_r_10_table);}
  if(mem_Glob_r_12_table)     		{clReleaseMemObject(mem_Glob_r_12_table);}
  if(mem_Glob_r_epsr_table)     	{clReleaseMemObject(mem_Glob_r_epsr_table);}
  if(mem_Glob_desolv_table)     	{clReleaseMemObject(mem_Glob_desolv_table);}
  if(mem_Glob_q1q2)     		{clReleaseMemObject(mem_Glob_q1q2);}
  if(mem_Glob_qasp_mul_absq)     	{clReleaseMemObject(mem_Glob_qasp_mul_absq);}
}

// Helper functions to display parameters returned by OpenCL queries
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name) {
   cl_ulong a;
   clGetDeviceInfo(device, param, sizeof(cl_ulong), &a, NULL);
   printf("%-40s = %lu\n", name, a);
}
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name) {
   cl_uint a;
   clGetDeviceInfo(device, param, sizeof(cl_uint), &a, NULL);
   printf("%-40s = %u\n", name, a);
}
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name) {
   cl_bool a;
   clGetDeviceInfo(device, param, sizeof(cl_bool), &a, NULL);
   printf("%-40s = %s\n", name, (a?"true":"false"));
}
static void device_info_string( cl_device_id device, cl_device_info param, const char* name) {
   char a[STRING_BUFFER_LEN]; 
   clGetDeviceInfo(device, param, STRING_BUFFER_LEN, &a, NULL);
   printf("%-40s = %s\n", name, a);
}

// Query and display OpenCL information on device and runtime environment
static void display_device_info( cl_device_id device ) {

   printf("Querying device for info:\n");
   printf("========================\n");
   device_info_string(device, CL_DEVICE_NAME, "CL_DEVICE_NAME");
   device_info_string(device, CL_DEVICE_VENDOR, "CL_DEVICE_VENDOR");
   device_info_uint(device, CL_DEVICE_VENDOR_ID, "CL_DEVICE_VENDOR_ID");
   device_info_string(device, CL_DEVICE_VERSION, "CL_DEVICE_VERSION");
   device_info_string(device, CL_DRIVER_VERSION, "CL_DRIVER_VERSION");
   device_info_uint(device, CL_DEVICE_ADDRESS_BITS, "CL_DEVICE_ADDRESS_BITS");
   device_info_bool(device, CL_DEVICE_AVAILABLE, "CL_DEVICE_AVAILABLE");
   device_info_bool(device, CL_DEVICE_ENDIAN_LITTLE, "CL_DEVICE_ENDIAN_LITTLE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_SIZE, "CL_DEVICE_GLOBAL_MEM_SIZE");
   device_info_bool(device, CL_DEVICE_IMAGE_SUPPORT, "CL_DEVICE_IMAGE_SUPPORT");
   device_info_ulong(device, CL_DEVICE_LOCAL_MEM_SIZE, "CL_DEVICE_LOCAL_MEM_SIZE");
   device_info_ulong(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, "CL_DEVICE_MAX_CLOCK_FREQUENCY");
   device_info_ulong(device, CL_DEVICE_MAX_COMPUTE_UNITS, "CL_DEVICE_MAX_COMPUTE_UNITS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_ARGS, "CL_DEVICE_MAX_CONSTANT_ARGS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE");
   device_info_uint(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, "CL_DEVICE_MEM_BASE_ADDR_ALIGN");
   device_info_uint(device, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, "CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE");

   {
      cl_command_queue_properties ccp;
      clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(cl_command_queue_properties), &ccp, NULL);
      printf("%-40s = %s\n", "Command queue out of order? ", ((ccp & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)?"true":"false"));
      printf("%-40s = %s\n", "Command queue profiling enabled? ", ((ccp & CL_QUEUE_PROFILING_ENABLE)?"true":"false"));
   }
}


