#include "performdocking.h"

#define STRINGIZE2(s) #s
#define STRINGIZE(s)	STRINGIZE2(s)
#define KRNL_BIN_FOLDER STRINGIZE(KRNL_LIB_DIRECTORY)
#define KRNL_SRC_FOLDER STRINGIZE(KRNL_SRC_DIRECTORY)
#define KRNL_SRC_COMMON STRINGIZE(KCMN_SRC_DIRECTORY)
#define KRNL0 	STRINGIZE(K0)
#define KRNL_GA STRINGIZE(K_GA)
#define KRNL_PC STRINGIZE(K_PC)
#define KRNL_IE STRINGIZE(K_IE)
#define KRNL_IA STRINGIZE(K_IA)
#define KRNL_LS STRINGIZE(K_LS)

//// --------------------------------
//// Device memory buffers
//// --------------------------------


/*								                  // Nr elements	// Nr bytes
cl_mem mem_atom_charges_const;		// float [MAX_NUM_OF_ATOMS];			// 90	 = 90	//360
cl_mem mem_atom_types_const;		// char  [MAX_NUM_OF_ATOMS];			// 90	 = 90	//360
cl_mem mem_intraE_contributors_const;	// char  [3*MAX_INTRAE_CONTRIBUTORS];		// 3*8128=28384 //24384	
cl_mem mem_VWpars_AC_const;		// float [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];	// 14*14 = 196  //784
cl_mem mem_VWpars_BD_const;		// float [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES]; // 14*14 = 196	//784
cl_mem mem_dspars_S_const;		// float [MAX_NUM_OF_ATYPES];			// 14    = 14   //56
cl_mem mem_dspars_V_const;		// float [MAX_NUM_OF_ATYPES];			// 14    = 14   //56
cl_mem mem_rotlist_const;		// int   [MAX_NUM_OF_ROTATIONS];		// 4096  = 4096 //16384
cl_mem mem_ref_coords_x_const;		// float [MAX_NUM_OF_ATOMS];			// 90    = 90   //360
cl_mem mem_ref_coords_y_const;		// float [MAX_NUM_OF_ATOMS];			// 90    = 90   //360
cl_mem mem_ref_coords_z_const;		// float [MAX_NUM_OF_ATOMS];			// 90	 = 90   //360
cl_mem mem_rotbonds_moving_vectors_const;// float [3*MAX_NUM_OF_ROTBONDS];		// 3*32  = 96   //384
cl_mem mem_rotbonds_unit_vectors_const;	// float [3*MAX_NUM_OF_ROTBONDS];		// 3*32  = 96   //384
cl_mem mem_ref_orientation_quats_const;	// float [4*MAX_NUM_OF_RUNS];			// 4*100 = 400  //1600
*/

//// --------------------------------
//// Docking
//// --------------------------------
int docking_with_fpga( const    Gridinfo*   mygrid,
	              /*const*/ float*      cpu_floatgrids,
                                Dockpars*   mypars,
		      const     Liganddata* myligand_init,
		      const     int*        argc,
		                char**      argv,
		                clock_t     clock_start_program)
/* The function performs the docking algorithm and generates the corresponding result files.
parameter mygrid:
		describes the grid
		filled with get_gridinfo()
parameter cpu_floatgrids:
		points to the memory region containing the grids
		filled with get_gridvalues_f()
parameter mypars:
		describes the docking parameters
		filled with get_commandpars()
parameter myligand_init:
		describes the ligands
		filled with get_liganddata()
parameters argc and argv:
		are the corresponding command line arguments parameter clock_start_program:
		contains the state of the clock tick counter at the beginning of the program
filled with clock() */
{
	// =======================================================================
	// Host Setup
	// =======================================================================

	const char* name_k_ga = KRNL_GA;
	const char* name_k_pc = KRNL_PC;
	const char* name_k_ie = KRNL_IE;
	const char* name_k_ia = KRNL_IA;
	const char* name_k_ls = KRNL_LS;

	const char* krnl_folder = KRNL_BIN_FOLDER;

	char path_k_ga[100];
	char path_k_pc[100];
	char path_k_ie[100];
	char path_k_ia[100];
	char path_k_ls[100];

	std::cout << "\n---------------------------------------------------------------------------------\n";
	std::cout << "Kernel binaries (VEO libraries)" << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n";
	strcpy(path_k_ga, krnl_folder); strcat(path_k_ga, "/"); strcat(path_k_ga, name_k_ga); strcat(path_k_ga, ".so"); std::cout << "path_k_ga: " << path_k_ga << std::endl;
	strcpy(path_k_pc, krnl_folder); strcat(path_k_pc, "/"); strcat(path_k_pc, name_k_pc); strcat(path_k_pc, ".so"); std::cout << "path_k_pc: " << path_k_pc << std::endl;
	strcpy(path_k_ie, krnl_folder); strcat(path_k_ie, "/"); strcat(path_k_ie, name_k_ie); strcat(path_k_ie, ".so"); std::cout << "path_k_ie: " << path_k_ie << std::endl;
	strcpy(path_k_ia, krnl_folder); strcat(path_k_ia, "/"); strcat(path_k_ia, name_k_ia); strcat(path_k_ia, ".so"); std::cout << "path_k_ia: " << path_k_ia << std::endl;
	strcpy(path_k_ls, krnl_folder); strcat(path_k_ls, "/"); strcat(path_k_ls, name_k_ls); strcat(path_k_ls, ".so"); std::cout << "path_k_ls: " << path_k_ls << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n" << std::endl;

	// VEO code
	wrapper_veo_api_version ();

	// Loading "ve_hello" on VE node 0
	struct veo_proc_handle *ve_process = wrapper_veo_proc_create(0);
	uint64_t kernel_ga_handle = wrapper_veo_load_library(ve_process, path_k_ga);
	uint64_t kernel_pc_handle = wrapper_veo_load_library(ve_process, path_k_pc);
	uint64_t kernel_ie_handle = wrapper_veo_load_library(ve_process, path_k_ie);
	uint64_t kernel_ia_handle = wrapper_veo_load_library(ve_process, path_k_ia);
	uint64_t kernel_ls_handle = wrapper_veo_load_library(ve_process, path_k_ls);
	struct veo_thr_ctxt *veo_thread_context = wrapper_veo_context_open(ve_process);

	// End of Host Setup
	// =======================================================================

	clock_t clock_start_docking;
	clock_t	clock_stop_docking;
	clock_t clock_stop_program_before_clustering;

	Liganddata myligand_reference;

 	// Allocating memory for floatgrids,
	size_t size_floatgrids_nbytes = sizeof(float) * (mygrid->num_of_atypes+2) *
					(mygrid->size_xyz[0]) * (mygrid->size_xyz[1]) * (mygrid->size_xyz[2]);

	size_t size_populations_nelems = mypars->num_of_runs * mypars->pop_size * ACTUAL_GENOTYPE_LENGTH;
	size_t size_populations_nbytes = size_populations_nelems * sizeof(float);
	
	size_t size_energies_nelems = mypars->num_of_runs * mypars->pop_size;
	size_t size_energies_nbytes = size_energies_nelems * sizeof(float);
	
	// Allocating and initializing CPU memory for initial population
	std::vector<float> cpu_init_populations (size_populations_nelems, 0.0f);

	// Allocating CPU memory for final population
	std::vector<float> cpu_final_populations (size_populations_nelems);

	// Allocating CPU memory for results
	std::vector<float> cpu_energies (size_energies_nelems);

	// Allocating CPU memory for resulting ligands
	std::vector<Ligandresult> cpu_result_ligands (mypars->num_of_runs);

	// Allocating memory in CPU for reference orientation angles
	std::vector<float> cpu_ref_ori_angles (mypars->num_of_runs*3);

	// Generating initial populations and random orientation angles of reference ligand
	// (ligand will be moved to origo and scaled as well)
	myligand_reference = *myligand_init;
	gen_initpop_and_reflig(mypars, cpu_init_populations.data(), cpu_ref_ori_angles.data(), &myligand_reference, mygrid);

	// Allocating memory in CPU for pseudorandom number generator seeds
	const unsigned int num_of_prng_blocks = 25;
	size_t size_prng_seeds_nelems = num_of_prng_blocks * mypars->num_of_runs;
	size_t size_prng_seeds_nbytes = size_prng_seeds_nelems * sizeof(unsigned int);
	std::vector<unsigned int> cpu_prng_seeds (size_prng_seeds_nelems);
	
	// Initializing seed generator
	genseed(time(NULL));	

	// Generating seeds (for each thread during GA)
	for (unsigned int i=0; i<size_prng_seeds_nelems; i++) {
		cpu_prng_seeds[i] = genseed(0u);
	}

	size_t size_evals_of_runs_nelems = mypars->num_of_runs;
	size_t size_evals_of_runs_nbytes = size_evals_of_runs_nelems * sizeof(int);

	// Allocating memory in CPU for evaluation counters
	std::vector<int> cpu_evals_of_runs (size_evals_of_runs_nelems, 0);
	
	// Allocating memory in CPU for generation counters
	std::vector<int> cpu_gens_of_runs (size_evals_of_runs_nelems, 0);

	// Preparing the constant data fields for the accelerator (calcenergy.cpp)
	// -----------------------------------------------------------------------------------------------------
	kernelconstant_static  KerConstStatic;
	if (prepare_conststatic_fields_for_fpga(&myligand_reference, mypars, cpu_ref_ori_angles.data(), &KerConstStatic) == 1)
		return 1;

	// Preparing parameter struct
	Dockparameters dockpars;
	dockpars.num_of_atoms  			= ((unsigned char)  myligand_reference.num_of_atoms);
	dockpars.num_of_atypes 			= ((unsigned char)  myligand_reference.num_of_atypes);
	dockpars.num_of_intraE_contributors 	= ((unsigned int)   myligand_reference.num_of_intraE_contributors);
	dockpars.gridsize_x    			= ((unsigned char)  mygrid->size_xyz[0]);
	dockpars.gridsize_y    			= ((unsigned char)  mygrid->size_xyz[1]);
	dockpars.gridsize_z    			= ((unsigned char)  mygrid->size_xyz[2]);
	dockpars.g1	       			= dockpars.gridsize_x ;
	dockpars.g2	       			= dockpars.gridsize_x * dockpars.gridsize_y;
	dockpars.g3	       			= dockpars.gridsize_x * dockpars.gridsize_y * dockpars.gridsize_z;
	dockpars.grid_spacing  			= ((float) mygrid->spacing);
	dockpars.rotbondlist_length 		= ((unsigned int) myligand_reference.num_of_rotcyc);
	dockpars.coeff_elec    			= ((float) mypars->coeffs.scaled_AD4_coeff_elec);
	dockpars.coeff_desolv  			= ((float) mypars->coeffs.AD4_coeff_desolv);
	dockpars.num_of_energy_evals 		= (unsigned int) mypars->num_of_energy_evals;
	dockpars.num_of_generations  		= (unsigned int) mypars->num_of_generations;
	dockpars.pop_size      			= (unsigned int) mypars->pop_size;
	dockpars.num_of_genes  			= (unsigned char)(myligand_reference.num_of_rotbonds + 6);
	dockpars.tournament_rate 		= (mypars->tournament_rate)/100;
	dockpars.crossover_rate  		= (mypars->crossover_rate)/100;
	dockpars.mutation_rate   		= (mypars->mutation_rate)/100;
	dockpars.abs_max_dang    		= mypars->abs_max_dang;
	dockpars.abs_max_dmov    		= mypars->abs_max_dmov;
	dockpars.lsearch_rate    		= mypars->lsearch_rate;
	dockpars.num_of_lsentities 		= (unsigned int) (mypars->lsearch_rate/100.0*mypars->pop_size + 0.5);
	dockpars.rho_lower_bound   		= mypars->rho_lower_bound;
	dockpars.base_dmov_mul_sqrt3 		= mypars->base_dmov_mul_sqrt3;
	dockpars.base_dang_mul_sqrt3 		= mypars->base_dang_mul_sqrt3;
	dockpars.cons_limit        		= (unsigned int) mypars->cons_limit;
	dockpars.max_num_of_iters  		= (unsigned int) mypars->max_num_of_iters;
	dockpars.qasp 				= mypars->qasp;
	dockpars.smooth 			= mypars->smooth;

	// These variables hold multiplications between kernel-constants
	// better calculate them here and then pass them to Krnl_GA
	const float two_absmaxdmov = 2.0 * dockpars.abs_max_dmov;
	const float two_absmaxdang = 2.0 * dockpars.abs_max_dang;

	// These variables hold multiplications between kernel-constants
	// better calculate them here and then pass them to Krnl_InterE and Krnl_InterE2
	const unsigned int mul_tmp2 = dockpars.num_of_atypes * dockpars.g3;
	const unsigned int mul_tmp3 = (dockpars.num_of_atypes + 1) * dockpars.g3;

	// -----------------------------------------------------------------------------------------------------
    // These commands will allocate memory on the FPGA. 

	// Krnl_GA buffers
/*	
	mem_dockpars_conformations_current_Initial (CL_MEM_READ_ONLY, size_populations_nbytes, cpu_init_populations.data());
	mem_dockpars_conformations_current_Final (CL_MEM_WRITE_ONLY, size_populations_nbytes,	cpu_final_populations.data());
	mem_dockpars_energies_current (CL_MEM_WRITE_ONLY, size_energies_nbytes, cpu_energies.data());
	mem_evals_performed (CL_MEM_WRITE_ONLY, size_evals_of_runs_nbytes, cpu_evals_of_runs.data());
	mem_gens_performed (CL_MEM_WRITE_ONLY, size_evals_of_runs_nbytes,	cpu_gens_of_runs.data());
*/
	uint64_t mem_dockpars_conformations_current_Initial;
	uint64_t mem_dockpars_conformations_current_Final;
	uint64_t mem_dockpars_energies_current;
	uint64_t mem_evals_performed;
	uint64_t mem_gens_performed;

	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_conformations_current_Initial, size_populations_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_conformations_current_Final, size_populations_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_energies_current, size_energies_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_evals_performed, size_evals_of_runs_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_gens_performed, size_evals_of_runs_nbytes);

	// Krnl_Conform buffers
/*	
	mem_KerConstStatic_rotlist_const	(CL_MEM_READ_ONLY, MAX_NUM_OF_ROTATIONS*sizeof(int), &KerConstStatic.rotlist_const[0]);
	mem_KerConstStatic_ref_coords_const	(CL_MEM_READ_ONLY, MAX_NUM_OF_ATOMS*sizeof(cl_float3), &KerConstStatic.ref_coords_const[0]);	
	mem_KerConstStatic_rotbonds_moving_vectors_const(CL_MEM_READ_ONLY, MAX_NUM_OF_ROTBONDS*sizeof(cl_float3), &KerConstStatic.rotbonds_moving_vectors_const[0]);
	mem_KerConstStatic_rotbonds_unit_vectors_const(CL_MEM_READ_ONLY, MAX_NUM_OF_ROTBONDS*sizeof(cl_float3),	&KerConstStatic.rotbonds_unit_vectors_const[0]);
	mem_KerConstStatic_ref_orientation_quats_const(CL_MEM_READ_ONLY, MAX_NUM_OF_RUNS*sizeof(cl_float4), &KerConstStatic.ref_orientation_quats_const[0]);
*/
	uint64_t mem_KerConstStatic_rotlist_const;
	uint64_t mem_KerConstStatic_ref_coords_const;
	uint64_t mem_KerConstStatic_rotbonds_moving_vectors_const;
	uint64_t mem_KerConstStatic_rotbonds_unit_vectors_const;
	uint64_t mem_KerConstStatic_ref_orientation_quats_const;

	size_t size_rotlist_nelems = MAX_NUM_OF_ROTATIONS;
	size_t size_rotlist_nbytes = size_rotlist_nelems * sizeof(int);

	size_t size_ref_coords_nelems = MAX_NUM_OF_ATOMS;
	size_t size_ref_coords_nbytes = size_ref_coords_nelems * 3 * sizeof(float);

	size_t size_rotbonds_moving_vectors_nelems = MAX_NUM_OF_ROTBONDS;
	size_t size_rotbonds_moving_vectors_nbytes = size_rotbonds_moving_vectors_nelems * 3 * sizeof(float);

	size_t size_rotbonds_unit_vectors_nelems = MAX_NUM_OF_ROTBONDS;
	size_t size_rotbonds_unit_vectors_nbytes = size_rotbonds_unit_vectors_nelems * 3 * sizeof(float);

	size_t size_ref_orientation_quats_nelems = MAX_NUM_OF_RUNS;
	size_t size_ref_orientation_quats_nbytes = size_ref_orientation_quats_nelems * 4 * sizeof(float);

	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_rotlist_const, size_rotlist_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_ref_coords_const, size_ref_coords_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_rotbonds_moving_vectors_const, size_rotbonds_moving_vectors_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_rotbonds_unit_vectors_const, size_rotbonds_unit_vectors_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_ref_orientation_quats_const, size_ref_orientation_quats_nbytes);

	// Krnl_InterE buffers
/*	
	mem_dockpars_fgrids (CL_MEM_READ_ONLY, size_floatgrids_nbytes, cpu_floatgrids);	
	mem_KerConstStatic_InterE_atom_charges_const(CL_MEM_READ_ONLY, MAX_NUM_OF_ATOMS*sizeof(float),	&KerConstStatic.atom_charges_const[0]);
	mem_KerConstStatic_InterE_atom_types_const(CL_MEM_READ_ONLY, MAX_NUM_OF_ATOMS*sizeof(char),	&KerConstStatic.atom_types_const[0]);
*/
	uint64_t mem_dockpars_fgrids;
	uint64_t mem_KerConstStatic_InterE_atom_charges_const;
	uint64_t mem_KerConstStatic_InterE_atom_types_const;

	size_t size_InterE_atom_charges_nelems = MAX_NUM_OF_ATOMS;
	size_t size_InterE_atom_charges_nbytes = size_InterE_atom_charges_nelems * sizeof(float);

	size_t size_InterE_atom_types_nelems = MAX_NUM_OF_ATOMS;
	size_t size_InterE_atom_types_nbytes = size_InterE_atom_types_nelems * sizeof(char);

	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_fgrids, size_floatgrids_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_InterE_atom_charges_const, size_InterE_atom_charges_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_InterE_atom_types_const, size_InterE_atom_types_nbytes);

	// Krnl_IntraE buffers
/*	
	mem_KerConstStatic_IntraE_atom_charges_const(CL_MEM_READ_ONLY, MAX_NUM_OF_ATOMS*sizeof(float),	&KerConstStatic.atom_charges_const[0]);
	mem_KerConstStatic_IntraE_atom_types_const(CL_MEM_READ_ONLY, MAX_NUM_OF_ATOMS*sizeof(char), &KerConstStatic.atom_types_const[0]);
	mem_KerConstStatic_intraE_contributors_const(CL_MEM_READ_ONLY, MAX_INTRAE_CONTRIBUTORS*sizeof(cl_char3), &KerConstStatic.intraE_contributors_const[0]);
	mem_KerConstStatic_reqm_const (CL_MEM_READ_ONLY, ATYPE_NUM*sizeof(float),	&KerConstStatic.reqm_const);	
	mem_KerConstStatic_reqm_hbond_const	(CL_MEM_READ_ONLY, ATYPE_NUM*sizeof(float),	&KerConstStatic.reqm_hbond_const);
	mem_KerConstStatic_atom1_types_reqm_const(CL_MEM_READ_ONLY, ATYPE_NUM*sizeof(unsigned int),	&KerConstStatic.atom1_types_reqm_const);
	mem_KerConstStatic_atom2_types_reqm_const(CL_MEM_READ_ONLY, ATYPE_NUM*sizeof(unsigned int),	&KerConstStatic.atom2_types_reqm_const);
	mem_KerConstStatic_VWpars_AC_const	(CL_MEM_READ_ONLY, MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES*sizeof(float),	&KerConstStatic.VWpars_AC_const[0]);
	mem_KerConstStatic_VWpars_BD_const	(CL_MEM_READ_ONLY, MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES*sizeof(float),	&KerConstStatic.VWpars_BD_const[0]);
	mem_KerConstStatic_dspars_S_const	(CL_MEM_READ_ONLY, MAX_NUM_OF_ATYPES*sizeof(float),			&KerConstStatic.dspars_S_const[0]);
	mem_KerConstStatic_dspars_V_const	(CL_MEM_READ_ONLY, MAX_NUM_OF_ATYPES*sizeof(float),			&KerConstStatic.dspars_V_const[0]);
*/
	uint64_t mem_KerConstStatic_IntraE_atom_charges_const;
	uint64_t mem_KerConstStatic_IntraE_atom_types_const;
	uint64_t mem_KerConstStatic_intraE_contributors_const;
	uint64_t mem_KerConstStatic_reqm_const;
	uint64_t mem_KerConstStatic_reqm_hbond_const;
	uint64_t mem_KerConstStatic_atom1_types_reqm_const;
	uint64_t mem_KerConstStatic_atom2_types_reqm_const;
	uint64_t mem_KerConstStatic_VWpars_AC_const;
	uint64_t mem_KerConstStatic_VWpars_BD_const;
	uint64_t mem_KerConstStatic_dspars_S_const;
	uint64_t mem_KerConstStatic_dspars_V_const;

	size_t size_IntraE_atom_charges_nelems = MAX_NUM_OF_ATOMS;
	size_t size_IntraE_atom_charges_nbytes = size_IntraE_atom_charges_nelems * sizeof(float);

	size_t size_IntraE_atom_types_nelems = MAX_NUM_OF_ATOMS;
	size_t size_IntraE_atom_types_nbytes = size_IntraE_atom_types_nelems * sizeof(char);

	size_t size_intraE_contributors_nelems = MAX_INTRAE_CONTRIBUTORS;
	size_t size_intraE_contributors_nbytes = size_intraE_contributors_nelems * 3 * sizeof(char);

	size_t size_reqm_nelems = ATYPE_NUM;
	size_t size_reqm_nbytes = size_reqm_nelems * sizeof(float);

	size_t size_reqm_hbond_nelems = ATYPE_NUM;
	size_t size_reqm_hbond_nbytes = size_reqm_hbond_nelems * sizeof(float);

	size_t size_atom1_types_reqm_nelems = ATYPE_NUM;
	size_t size_atom1_types_reqm_nbytes = size_atom1_types_reqm_nelems * sizeof(unsigned int);

	size_t size_atom2_types_reqm_nelems = ATYPE_NUM;
	size_t size_atom2_types_reqm_nbytes = size_atom2_types_reqm_nelems * sizeof(unsigned int);

	size_t size_VWpars_AC_nelems = MAX_NUM_OF_ATYPES * MAX_NUM_OF_ATYPES;
	size_t size_VWpars_AC_nbytes = size_VWpars_AC_nelems * sizeof(float);

	size_t size_VWpars_BD_nelems = MAX_NUM_OF_ATYPES * MAX_NUM_OF_ATYPES;
	size_t size_VWpars_BD_nbytes = size_VWpars_BD_nelems * sizeof(float);

	size_t size_dspars_S_nelems = MAX_NUM_OF_ATYPES;
	size_t size_dspars_S_nbytes = size_dspars_S_nelems * sizeof(float);

	size_t size_dspars_V_nelems = MAX_NUM_OF_ATYPES;
	size_t size_dspars_V_nbytes = size_dspars_V_nelems * sizeof(float);

	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_IntraE_atom_charges_const, size_IntraE_atom_charges_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_IntraE_atom_types_const, size_IntraE_atom_types_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_intraE_contributors_const, size_intraE_contributors_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_reqm_const, size_reqm_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_reqm_hbond_const, size_reqm_hbond_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_atom1_types_reqm_const, size_atom1_types_reqm_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_atom2_types_reqm_const, size_atom2_types_reqm_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_VWpars_AC_const, size_VWpars_AC_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_VWpars_BD_const, size_VWpars_BD_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_dspars_S_const, size_dspars_S_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_KerConstStatic_dspars_V_const, size_dspars_V_nbytes);

	// FIXME: should be commented out for the final version ?
	// Printing sizes
	std::cout << "\n---------------------------------------------------------------------------------\n";
	std::cout << std::left << std::setw(40) << "Memory sizes" << std::right << std::setw(31) << "Bytes" << std::right << std::setw(10) << "KB" << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n";

	std::cout << std::left << std::setw(40) << "size_floatgrids_nbytes" << std::right << std::setw(31) << size_floatgrids_nbytes << std::right << std::setw(10) << sizeKB(size_floatgrids_nbytes) << std::endl;
	std::cout << std::left << std::setw(40) << "size_populations_nbytes" << std::right << std::setw(31) << size_populations_nbytes << std::right << std::setw(10) << sizeKB(size_populations_nbytes) << std::endl;
	std::cout << std::left << std::setw(40) << "size_energies_nbytes" << std::right << std::setw(31) << size_energies_nbytes << std::right << std::setw(10) << sizeKB(size_energies_nbytes) << std::endl;
	std::cout << std::left << std::setw(40) << "size_evals_of_runs_nbytes" << std::right << std::setw(31) << size_evals_of_runs_nbytes << std::right << std::setw(10) << sizeKB(size_evals_of_runs_nbytes) << std::endl;
	std::cout << std::left << std::setw(40) << "size_prng_seeds_nbytes" << std::right << std::setw(31) << size_prng_seeds_nbytes << std::right << std::setw(10) << sizeKB(size_prng_seeds_nbytes) << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n" << std::endl;


	// -----------------------------------------------------------------------------------------------------
/*
	//Separate Read/write Buffer vector is needed to migrate data between host/device
	std::vector<cl::Memory> inBufVec, outBufVec;

	// Krnl_GA
	inBufVec.push_back(mem_dockpars_conformations_current_Initial);
	// Krnl_Conform
	inBufVec.push_back(mem_KerConstStatic_rotlist_const);
	inBufVec.push_back(mem_KerConstStatic_ref_coords_const);
	inBufVec.push_back(mem_KerConstStatic_rotbonds_moving_vectors_const);
	inBufVec.push_back(mem_KerConstStatic_rotbonds_unit_vectors_const);
	inBufVec.push_back(mem_KerConstStatic_ref_orientation_quats_const);
	// Krnl_InterE
	inBufVec.push_back(mem_dockpars_fgrids);
	inBufVec.push_back(mem_KerConstStatic_InterE_atom_charges_const);
	inBufVec.push_back(mem_KerConstStatic_InterE_atom_types_const);
	// Krnl_IntraE
	inBufVec.push_back(mem_KerConstStatic_IntraE_atom_charges_const);
	inBufVec.push_back(mem_KerConstStatic_IntraE_atom_types_const);
	inBufVec.push_back(mem_KerConstStatic_intraE_contributors_const);
	inBufVec.push_back(mem_KerConstStatic_reqm_const);
	inBufVec.push_back(mem_KerConstStatic_reqm_hbond_const);
	inBufVec.push_back(mem_KerConstStatic_atom1_types_reqm_const);
	inBufVec.push_back(mem_KerConstStatic_atom2_types_reqm_const);
	inBufVec.push_back(mem_KerConstStatic_VWpars_AC_const);
	inBufVec.push_back(mem_KerConstStatic_VWpars_BD_const);
	inBufVec.push_back(mem_KerConstStatic_dspars_S_const);
	inBufVec.push_back(mem_KerConstStatic_dspars_V_const);

	// Krnl_GA
    outBufVec.push_back(mem_dockpars_conformations_current_Final);
	outBufVec.push_back(mem_dockpars_energies_current);
	outBufVec.push_back(mem_evals_performed);
	outBufVec.push_back(mem_gens_performed);
*/
	// -----------------------------------------------------------------------------------------------------
	// These commands will load CPU-sources vectors from the host
   	// application and into cl::Buffer objects. 
	// The data will be be transferred from system memory 
	// over PCIe to the FPGA on-board DDR memory.
/*	
	command_queue_ga.enqueueMigrateMemObjects(inBufVec,0); // 2nd arg 0 means from host
*/	



	// -----------------------------------------------------------------------------------------------------

	clock_start_docking = clock();

	int narg;
/*	
	#ifdef ENABLE_KRNL_GA
*/	
	narg = 0;
	kernel_ga.setArg(narg++, mem_dockpars_conformations_current_Initial);
	kernel_ga.setArg(narg++, mem_dockpars_conformations_current_Final);
	kernel_ga.setArg(narg++, mem_dockpars_energies_current);
	kernel_ga.setArg(narg++, mem_evals_performed);
	kernel_ga.setArg(narg++, mem_gens_performed);
	kernel_ga.setArg(narg++, dockpars.pop_size);
	kernel_ga.setArg(narg++, dockpars.num_of_energy_evals);
	kernel_ga.setArg(narg++, dockpars.num_of_generations);
	kernel_ga.setArg(narg++, dockpars.tournament_rate);
	kernel_ga.setArg(narg++, dockpars.mutation_rate);
	kernel_ga.setArg(narg++, dockpars.abs_max_dmov);
	kernel_ga.setArg(narg++, dockpars.abs_max_dang);
	kernel_ga.setArg(narg++, two_absmaxdmov);
	kernel_ga.setArg(narg++, two_absmaxdang);
	kernel_ga.setArg(narg++, dockpars.crossover_rate);
	kernel_ga.setArg(narg++, dockpars.num_of_lsentities);
	kernel_ga.setArg(narg++, dockpars.num_of_genes);
	// Other kernel args are configured at every docking run
/*	
	#endif
	#ifdef ENABLE_KRNL_CONFORM
*/	
	narg = 0;
	kernel_conform.setArg(narg++, mem_KerConstStatic_rotlist_const);
	kernel_conform.setArg(narg++, mem_KerConstStatic_ref_coords_const);
	kernel_conform.setArg(narg++, mem_KerConstStatic_rotbonds_moving_vectors_const);
	kernel_conform.setArg(narg++, mem_KerConstStatic_rotbonds_unit_vectors_const);
	kernel_conform.setArg(narg++, dockpars.rotbondlist_length);
	kernel_conform.setArg(narg++, dockpars.num_of_atoms);
	kernel_conform.setArg(narg++, dockpars.num_of_genes);
	kernel_conform.setArg(narg++, mem_KerConstStatic_ref_orientation_quats_const);
	// Other kernel args are configured at every docking run
/*	
	#endif
*/
	unsigned char gridsizex_minus1 = dockpars.gridsize_x - 1;
	unsigned char gridsizey_minus1 = dockpars.gridsize_y - 1;
	unsigned char gridsizez_minus1 = dockpars.gridsize_z - 1;
	float fgridsizex_minus1 = (float) gridsizex_minus1;
	float fgridsizey_minus1 = (float) gridsizey_minus1;
	float fgridsizez_minus1 = (float) gridsizez_minus1;
/*
	#ifdef ENABLE_KRNL_INTERE
*/	
	narg = 0;
	kernel_intere.setArg(narg++, mem_dockpars_fgrids);
	kernel_intere.setArg(narg++, mem_KerConstStatic_InterE_atom_charges_const);
	kernel_intere.setArg(narg++, mem_KerConstStatic_InterE_atom_types_const);
	kernel_intere.setArg(narg++, dockpars.g1);
	kernel_intere.setArg(narg++, dockpars.g2);
	kernel_intere.setArg(narg++, dockpars.g3);
	kernel_intere.setArg(narg++, dockpars.num_of_atoms);
	kernel_intere.setArg(narg++, fgridsizex_minus1);
	kernel_intere.setArg(narg++, fgridsizey_minus1);
	kernel_intere.setArg(narg++, fgridsizez_minus1);
	kernel_intere.setArg(narg++, mul_tmp2);
	kernel_intere.setArg(narg++, mul_tmp3);
/*	
	#endif
	#ifdef ENABLE_KRNL_INTRAE
*/	
	narg = 0;
	kernel_intrae.setArg(narg++, mem_KerConstStatic_IntraE_atom_charges_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_IntraE_atom_types_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_intraE_contributors_const);
	kernel_intrae.setArg(narg++, dockpars.smooth);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_reqm_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_reqm_hbond_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_atom1_types_reqm_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_atom2_types_reqm_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_VWpars_AC_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_VWpars_BD_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_dspars_S_const);
	kernel_intrae.setArg(narg++, mem_KerConstStatic_dspars_V_const);
	kernel_intrae.setArg(narg++, dockpars.num_of_atoms);
	kernel_intrae.setArg(narg++, dockpars.num_of_intraE_contributors);
	kernel_intrae.setArg(narg++, dockpars.grid_spacing);
	kernel_intrae.setArg(narg++, dockpars.num_of_atypes);
	kernel_intrae.setArg(narg++, dockpars.coeff_elec);
	kernel_intrae.setArg(narg++, dockpars.qasp);
	kernel_intrae.setArg(narg++, dockpars.coeff_desolv);
/*	
	#endif
*/
	unsigned short Host_max_num_of_iters = (unsigned short)dockpars.max_num_of_iters;
	unsigned char  Host_cons_limit       = (unsigned char) dockpars.cons_limit;
/*
	#ifdef ENABLE_KRNL_LS
*/	
	narg = 0;
	kernel_ls.setArg(narg++, Host_max_num_of_iters);
	kernel_ls.setArg(narg++, dockpars.rho_lower_bound);
	kernel_ls.setArg(narg++, dockpars.base_dmov_mul_sqrt3);
	kernel_ls.setArg(narg++, dockpars.num_of_genes);
	kernel_ls.setArg(narg++, dockpars.base_dang_mul_sqrt3);
	kernel_ls.setArg(narg++, Host_cons_limit);
/*	
	#endif
*/
	printf("Docking runs to be executed: %lu\n", mypars->num_of_runs); 
	printf("Execution run: ");

	for (unsigned int run_cnt = 0; run_cnt < mypars->num_of_runs; run_cnt++) {
		printf(" %u", run_cnt+1); 
		fflush(stdout);
/*
		#ifdef ENABLE_KRNL_GA
*/		
		unsigned short ushort_run_cnt  = (unsigned ushort) run_cnt;
		unsigned int   Host_Offset_Pop = run_cnt * dockpars.pop_size * ACTUAL_GENOTYPE_LENGTH;
		unsigned int   Host_Offset_Ene = run_cnt * dockpars.pop_size;
		kernel_ga.setArg(17, ushort_run_cnt);
		kernel_ga.setArg(18, Host_Offset_Pop);
		kernel_ga.setArg(19, Host_Offset_Ene);
/*		
		#endif
		#ifdef ENABLE_KRNL_CONFORM
*/		
		kernel_conform.setArg(8, ushort_run_cnt);
/*		
		#endif
		#ifdef ENABLE_KRNL_PRNG_BT_USHORT_FLOAT
*/		
		kernel_prng_bt_ushort_float.setArg(0, cpu_prng_seeds[num_of_prng_blocks * run_cnt]);
		kernel_prng_bt_ushort_float.setArg(1, cpu_prng_seeds[num_of_prng_blocks * run_cnt + 1]);
/*		
		#endif
		#ifdef ENABLE_KRNL_GA
*/		
		command_queue_ga.enqueueTask(kernel_ga);
/*		
		#endif
		#ifdef ENABLE_KRNL_CONFORM
*/		
		command_queue_conform.enqueueTask(kernel_conform);
/*		
		#endif
		#ifdef ENABLE_KRNL_INTERE
*/		
		command_queue_intere.enqueueTask(kernel_intere);
/*		
		#endif
		#ifdef ENABLE_KRNL_INTRAE
*/		
		command_queue_intrae.enqueueTask(kernel_intrae);
/*		
		#endif
		#ifdef ENABLE_KRNL_PRNG_BT_USHORT_FLOAT
*/		
		command_queue_prng_bt_ushort_float.enqueueTask(kernel_prng_bt_ushort_float);
/*		
		#endif
		#ifdef ENABLE_KRNL_LS
*/		
		command_queue_ls.enqueueTask(kernel_ls);
/*		
		#endif
*/
		clock_stop_docking = clock();
	} // End of for (unsigned int run_cnt = 0; run_cnt < mypars->num_of_runs; run_cnt++)

	printf("\n");
	fflush(stdout);

	// -----------------------------------------------------------------------------------------------------
    	// These commands will load CPU-sources vectors from the host
   	// application and into cl::Buffer objects. 
	// The data will be be transferred from system memory 
	// over PCIe to the FPGA on-board DDR memory.
/*	
	command_queue_ga.enqueueMigrateMemObjects(outBufVec,CL_MIGRATE_MEM_OBJECT_HOST);
*/
	// -----------------------------------------------------------------------------------------------------

	for (unsigned int run_cnt = 0; run_cnt < mypars->num_of_runs; run_cnt++) {

		arrange_result(cpu_final_populations.data() + run_cnt*mypars->pop_size*ACTUAL_GENOTYPE_LENGTH, 
			       cpu_energies.data()          + run_cnt*mypars->pop_size, 
			       mypars->pop_size);

		/*printf("cpu_evals_of_runs[%u]: %u\n", run_cnt, cpu_evals_of_runs[run_cnt]);*/

		make_resfiles(cpu_final_populations.data() + run_cnt*mypars->pop_size*ACTUAL_GENOTYPE_LENGTH, 
			      cpu_energies.data()          + run_cnt*mypars->pop_size, 
			      &myligand_reference,
			      myligand_init, 
			      mypars, 
   			      cpu_evals_of_runs[run_cnt], 
			      cpu_gens_of_runs[run_cnt], /*generation_cnt, */
			      mygrid, 
			      cpu_floatgrids,
			      cpu_ref_ori_angles.data() + 3*run_cnt, 
			      argc, 
			      argv, 
			      0,
			      run_cnt, 
                              &(cpu_result_ligands [run_cnt]));
	} // End of for (unsigned int run_cnt = 0; run_cnt < mypars->num_of_runs; run_cnt++) 

	/*
	for (int cnt_pop=0;cnt_pop<size_populations/sizeof(float);cnt_pop++)
		printf("total_num_pop: %u, cpu_final_populations[%u]: %f\n",(unsigned int)(size_populations/sizeof(float)),cnt_pop,cpu_final_populations[cnt_pop]);

	for (int cnt_pop=0;cnt_pop<size_energies/sizeof(float);cnt_pop++)
		printf("total_num_energies: %u, cpu_energies[%u]: %f\n",    (unsigned int)(size_energies/sizeof(float)),cnt_pop,cpu_energies[cnt_pop]);
	*/

	clock_stop_program_before_clustering = clock();


	clusanal_gendlg(cpu_result_ligands.data(), 
			mypars->num_of_runs,
			myligand_init, mypars,
   		        mygrid, 
			argc,
			argv,
			ELAPSEDSECS(clock_stop_docking, clock_start_docking)/mypars->num_of_runs,
			ELAPSEDSECS(clock_stop_program_before_clustering, clock_start_program));

	clock_stop_docking = clock();

	return 0;
}
