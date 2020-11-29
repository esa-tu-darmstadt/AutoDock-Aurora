#include "performdocking.h"

#define STRINGIZE2(s) #s
#define STRINGIZE(s)	STRINGIZE2(s)
#define KRNL_BIN_FOLDER STRINGIZE(KRNL_LIB_DIRECTORY)
#define KRNL_SRC_FOLDER STRINGIZE(KRNL_SRC_DIRECTORY)
#define KRNL_SRC_COMMON STRINGIZE(KCMN_SRC_DIRECTORY)
#define KRNL0 	STRINGIZE(K0)
#define KRNL_GA STRINGIZE(K_GA)

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
int docking_with_aurora(
	const	Gridinfo*	mygrid,
	/*const*/ float*	cpu_floatgrids,
			Dockpars*	mypars,
	const	Liganddata*	myligand_init,
	const 	Liganddata* myxrayligand,
	const	int*		argc,
			char**		argv,
			clock_t		clock_start_program
)
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
parameter myxrayligand:
		describes the xray ligand
		filled with get_xrayliganddata()
parameters argc and argv:
		are the corresponding command line arguments parameter clock_start_program:
		contains the state of the clock tick counter at the beginning of the program
filled with clock() */
{
	// =======================================================================
	// Host Setup
	// =======================================================================

	const char* name_k_ga = KRNL_GA;
	const char* krnl_folder = KRNL_BIN_FOLDER;
	char path_k_ga[100];

	std::cout << "\n---------------------------------------------------------------------------------\n";
	std::cout << "Kernel binaries (VEO libraries)" << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n";
	strcpy(path_k_ga, krnl_folder); strcat(path_k_ga, "/"); strcat(path_k_ga, name_k_ga); strcat(path_k_ga, ".so"); std::cout << "path_k_ga: " << path_k_ga << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n" << std::endl;

	// VEO code
	wrapper_veo_api_version ();

	// Loading "ve_hello" on VE node 0
	struct veo_proc_handle *ve_process = wrapper_veo_proc_create(0);
	uint64_t kernel_ga_handle = wrapper_veo_load_library(ve_process, path_k_ga);
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

	// TODO: check passing of random numbers
	// TODO: use parameter instead of hardcoding num_genes
	// Allocating memory in CPU for pseudorandom number generator seeds
	size_t size_prng_seeds_nelems = mypars->num_of_runs * mypars->pop_size;
	size_t size_prng_seeds_nbytes = size_prng_seeds_nelems * sizeof(unsigned int);
	std::vector<unsigned int> cpu_prng_seeds (size_prng_seeds_nelems);
	
	// Initializing seed generator
	genseed(time(NULL));	

	// Generating seeds (for each thread during GA)
	for (unsigned int i = 0; i < size_prng_seeds_nelems; i++) {
		cpu_prng_seeds[i] = genseed(0u);
	}

	size_t size_evals_of_runs_nelems = mypars->num_of_runs;
	size_t size_evals_of_runs_nbytes = size_evals_of_runs_nelems * sizeof(int);

	// Allocating memory in CPU for evaluation counters
	std::vector<int> cpu_evals_of_runs (size_evals_of_runs_nelems, 0);
	
	// Allocating memory in CPU for generation counters
	std::vector<int> cpu_gens_of_runs (size_evals_of_runs_nelems, 0);

	// -----------------------------------------------------------------------------------------------------
	// Preparing the constant data fields for the accelerator (calcenergy.cpp)
	// -----------------------------------------------------------------------------------------------------

	kernelconstant_static  KerConstStatic;
	if (prepare_conststatic_fields_for_fpga(&myligand_reference, mypars, cpu_ref_ori_angles.data(), &KerConstStatic) == 1)
		return 1;

	// Preparing parameter struct
	Dockparameters dockpars;
	dockpars.num_of_atoms  			= ((unsigned char)  myligand_reference.num_of_atoms);
	dockpars.num_of_atypes 			= ((unsigned int)   myligand_reference.num_of_atypes);
	dockpars.num_of_intraE_contributors 	= ((unsigned int)   myligand_reference.num_of_intraE_contributors);
	dockpars.gridsize_x    			= ((unsigned char)  mygrid->size_xyz[0]);
	dockpars.gridsize_y    			= ((unsigned char)  mygrid->size_xyz[1]);
	dockpars.gridsize_z    			= ((unsigned char)  mygrid->size_xyz[2]);
	dockpars.g1	       				= dockpars.gridsize_x ;
	dockpars.g2	       				= dockpars.gridsize_x * dockpars.gridsize_y;
	dockpars.g3	       				= dockpars.gridsize_x * dockpars.gridsize_y * dockpars.gridsize_z;
	dockpars.grid_spacing  			= ((float) mygrid->spacing);
	dockpars.rotbondlist_length 	= ((unsigned int) myligand_reference.num_of_rotcyc);
	dockpars.coeff_elec    			= ((float) mypars->coeffs.scaled_AD4_coeff_elec);
	dockpars.coeff_desolv  			= ((float) mypars->coeffs.AD4_coeff_desolv);
	dockpars.num_of_energy_evals 	= (unsigned int) mypars->num_of_energy_evals;
	dockpars.num_of_generations  	= (unsigned int) mypars->num_of_generations;
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
	dockpars.base_dmov_mul_sqrt3 	= mypars->base_dmov_mul_sqrt3;
	dockpars.base_dang_mul_sqrt3 	= mypars->base_dang_mul_sqrt3;
	dockpars.cons_limit        		= (unsigned int) mypars->cons_limit;
	dockpars.max_num_of_iters  		= (unsigned int) mypars->max_num_of_iters;
	dockpars.qasp 					= mypars->qasp;
	dockpars.smooth 				= mypars->smooth;

	// These variables hold multiplications between kernel-constants
	// better calculate them here and then pass them to Krnl_GA
	const float two_absmaxdmov = 2.0 * dockpars.abs_max_dmov;
	const float two_absmaxdang = 2.0 * dockpars.abs_max_dang;

	// These variables hold multiplications between kernel-constants
	// better calculate them here and then pass them to Krnl_InterE and Krnl_InterE2
	const unsigned int mul_tmp2 = dockpars.num_of_atypes * dockpars.g3;
	const unsigned int mul_tmp3 = (dockpars.num_of_atypes + 1) * dockpars.g3;

	// -----------------------------------------------------------------------------------------------------
	// Printing sizes
	// -----------------------------------------------------------------------------------------------------

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
	// Defining kernel buffers
	// -----------------------------------------------------------------------------------------------------

	uint64_t mem_dockpars_conformations_current_Initial;
	uint64_t mem_dockpars_conformations_current_Final;
	uint64_t mem_dockpars_energies_current;
	uint64_t mem_evals_performed;
	uint64_t mem_gens_performed;
	uint64_t mem_prng_states;

	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_conformations_current_Initial, size_populations_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_conformations_current_Final, size_populations_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_energies_current, size_energies_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_evals_performed, size_evals_of_runs_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_gens_performed, size_evals_of_runs_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_prng_states, size_prng_seeds_nbytes);

	wrapper_veo_write_mem (ve_process, mem_dockpars_conformations_current_Initial, cpu_init_populations.data(), size_populations_nbytes);
	wrapper_veo_write_mem (ve_process, mem_prng_states, cpu_prng_seeds.data(), size_prng_seeds_nbytes);

	// Pose Calculation buffers
	uint64_t mem_pc_rotlist_const;
	uint64_t mem_pc_ref_coords_x_const;	// TODO: fusion x,y,z into a single one
	uint64_t mem_pc_ref_coords_y_const;
	uint64_t mem_pc_ref_coords_z_const;
	uint64_t mem_pc_rotbonds_moving_vectors_const;
	uint64_t mem_pc_rotbonds_unit_vectors_const;
	uint64_t mem_pc_ref_orientation_quats_const;

	size_t size_rotlist_nelems = MAX_NUM_OF_ROTATIONS;
	size_t size_rotlist_nbytes = size_rotlist_nelems * sizeof(int);

	size_t size_ref_coords_nelems = MAX_NUM_OF_ATOMS;
	size_t size_ref_coords_nbytes = size_ref_coords_nelems * sizeof(float);

	size_t size_rotbonds_moving_vectors_nelems = MAX_NUM_OF_ROTBONDS;
	size_t size_rotbonds_moving_vectors_nbytes = size_rotbonds_moving_vectors_nelems * 3 * sizeof(float);

	size_t size_rotbonds_unit_vectors_nelems = MAX_NUM_OF_ROTBONDS;
	size_t size_rotbonds_unit_vectors_nbytes = size_rotbonds_unit_vectors_nelems * 3 * sizeof(float);

	size_t size_ref_orientation_quats_nelems = MAX_NUM_OF_RUNS;
	size_t size_ref_orientation_quats_nbytes = size_ref_orientation_quats_nelems * 4 * sizeof(float);

	wrapper_veo_alloc_mem (ve_process, &mem_pc_rotlist_const, size_rotlist_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_pc_ref_coords_x_const, size_ref_coords_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_pc_ref_coords_y_const, size_ref_coords_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_pc_ref_coords_z_const, size_ref_coords_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_pc_rotbonds_moving_vectors_const, size_rotbonds_moving_vectors_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_pc_rotbonds_unit_vectors_const, size_rotbonds_unit_vectors_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_pc_ref_orientation_quats_const, size_ref_orientation_quats_nbytes);

	wrapper_veo_write_mem (ve_process, mem_pc_rotlist_const, &KerConstStatic.rotlist_const[0], size_rotlist_nbytes);
	wrapper_veo_write_mem (ve_process, mem_pc_ref_coords_x_const, &KerConstStatic.ref_coords_x_const[0], size_ref_coords_nbytes);
	wrapper_veo_write_mem (ve_process, mem_pc_ref_coords_y_const, &KerConstStatic.ref_coords_y_const[0], size_ref_coords_nbytes);
	wrapper_veo_write_mem (ve_process, mem_pc_ref_coords_z_const, &KerConstStatic.ref_coords_z_const[0], size_ref_coords_nbytes);
	wrapper_veo_write_mem (ve_process, mem_pc_rotbonds_moving_vectors_const, &KerConstStatic.rotbonds_moving_vectors_const[0], size_rotbonds_moving_vectors_nbytes);
	wrapper_veo_write_mem (ve_process, mem_pc_rotbonds_unit_vectors_const, &KerConstStatic.rotbonds_unit_vectors_const[0], size_rotbonds_unit_vectors_nbytes);
	wrapper_veo_write_mem (ve_process, mem_pc_ref_orientation_quats_const, &KerConstStatic.ref_orientation_quats_const[0], size_ref_orientation_quats_nbytes);

	// IE buffers
	uint64_t mem_dockpars_fgrids;
	uint64_t mem_ia_ie_atom_charges_const;
	uint64_t mem_ia_ie_atom_types_const;

	size_t size_InterE_atom_charges_nelems = MAX_NUM_OF_ATOMS;
	size_t size_InterE_atom_charges_nbytes = size_InterE_atom_charges_nelems * sizeof(float);

	size_t size_InterE_atom_types_nelems = MAX_NUM_OF_ATOMS;
	size_t size_InterE_atom_types_nbytes = size_InterE_atom_types_nelems * sizeof(char);

	wrapper_veo_alloc_mem (ve_process, &mem_dockpars_fgrids, size_floatgrids_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_ie_atom_charges_const, size_InterE_atom_charges_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_ie_atom_types_const, size_InterE_atom_types_nbytes);

	wrapper_veo_write_mem (ve_process, mem_dockpars_fgrids, cpu_floatgrids, size_floatgrids_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_ie_atom_charges_const, &KerConstStatic.atom_charges_const[0], size_InterE_atom_charges_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_ie_atom_types_const, &KerConstStatic.atom_types_const[0], size_InterE_atom_types_nbytes);

	// IA buffers
	uint64_t mem_ia_contributors_const;
	uint64_t mem_ia_reqm_const;
	uint64_t mem_ia_reqm_hbond_const;
	uint64_t mem_ia_atom1_types_reqm_const;
	uint64_t mem_ia_atom2_types_reqm_const;
	uint64_t mem_ia_VWpars_AC_const;
	uint64_t mem_ia_VWpars_BD_const;
	uint64_t mem_ia_dspars_S_const;
	uint64_t mem_ia_dspars_V_const;

	size_t size_IntraE_atom_charges_nelems = MAX_NUM_OF_ATOMS;
	size_t size_IntraE_atom_charges_nbytes = size_IntraE_atom_charges_nelems * sizeof(float);

	size_t size_IntraE_atom_types_nelems = MAX_NUM_OF_ATOMS;
	size_t size_IntraE_atom_types_nbytes = size_IntraE_atom_types_nelems * sizeof(int);

	size_t size_intraE_contributors_nelems = MAX_INTRAE_CONTRIBUTORS;
	size_t size_intraE_contributors_nbytes = size_intraE_contributors_nelems * 3 * sizeof(int);

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

	wrapper_veo_alloc_mem (ve_process, &mem_ia_contributors_const, size_intraE_contributors_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_reqm_const, size_reqm_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_reqm_hbond_const, size_reqm_hbond_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_atom1_types_reqm_const, size_atom1_types_reqm_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_atom2_types_reqm_const, size_atom2_types_reqm_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_VWpars_AC_const, size_VWpars_AC_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_VWpars_BD_const, size_VWpars_BD_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_dspars_S_const, size_dspars_S_nbytes);
	wrapper_veo_alloc_mem (ve_process, &mem_ia_dspars_V_const, size_dspars_V_nbytes);

	wrapper_veo_write_mem (ve_process, mem_ia_contributors_const, &KerConstStatic.intraE_contributors_const[0], size_intraE_contributors_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_reqm_const, &KerConstStatic.reqm_const, size_reqm_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_reqm_hbond_const, &KerConstStatic.reqm_hbond_const, size_reqm_hbond_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_atom1_types_reqm_const, &KerConstStatic.atom1_types_reqm_const, size_atom1_types_reqm_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_atom2_types_reqm_const, &KerConstStatic.atom2_types_reqm_const, size_atom2_types_reqm_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_VWpars_AC_const, &KerConstStatic.VWpars_AC_const[0], size_VWpars_AC_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_VWpars_BD_const, &KerConstStatic.VWpars_BD_const[0], size_VWpars_BD_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_dspars_S_const, &KerConstStatic.dspars_S_const[0], size_dspars_S_nbytes);
	wrapper_veo_write_mem (ve_process, mem_ia_dspars_V_const, &KerConstStatic.dspars_V_const[0], size_dspars_V_nbytes);

	// -----------------------------------------------------------------------------------------------------
	// Defining kernel arguments
	// -----------------------------------------------------------------------------------------------------

	// Defining extra variables

	// IE
	unsigned char gridsizex_minus1 = dockpars.gridsize_x - 1;
	unsigned char gridsizey_minus1 = dockpars.gridsize_y - 1;
	unsigned char gridsizez_minus1 = dockpars.gridsize_z - 1;
	float fgridsizex_minus1 = (float) gridsizex_minus1;
	float fgridsizey_minus1 = (float) gridsizey_minus1;
	float fgridsizez_minus1 = (float) gridsizez_minus1;

	// LS
	unsigned short Host_max_num_of_iters = (unsigned short)dockpars.max_num_of_iters;
	unsigned char  Host_cons_limit       = (unsigned char) dockpars.cons_limit;

	// Creating a VEO arguments object
	int narg = 0;
	struct veo_args *kernel_ga_arg_ptr = wrapper_veo_args_alloc ();

	// GA
	wrapper_veo_args_set_u64   (kernel_ga_arg_ptr, narg++, mem_dockpars_conformations_current_Initial);
	wrapper_veo_args_set_u64   (kernel_ga_arg_ptr, narg++, mem_dockpars_conformations_current_Final);
	wrapper_veo_args_set_u64   (kernel_ga_arg_ptr, narg++, mem_dockpars_energies_current);
	wrapper_veo_args_set_u64   (kernel_ga_arg_ptr, narg++, mem_evals_performed);
	wrapper_veo_args_set_u64   (kernel_ga_arg_ptr, narg++, mem_gens_performed);
	wrapper_veo_args_set_u64   (kernel_ga_arg_ptr, narg++, mem_prng_states);
	wrapper_veo_args_set_i32   (kernel_ga_arg_ptr, narg++, dockpars.pop_size);
	wrapper_veo_args_set_u32   (kernel_ga_arg_ptr, narg++, dockpars.num_of_energy_evals);
	wrapper_veo_args_set_u32   (kernel_ga_arg_ptr, narg++, dockpars.num_of_generations);
	wrapper_veo_args_set_float (kernel_ga_arg_ptr, narg++, dockpars.tournament_rate);
	wrapper_veo_args_set_float (kernel_ga_arg_ptr, narg++, dockpars.mutation_rate);
	wrapper_veo_args_set_float (kernel_ga_arg_ptr, narg++, dockpars.abs_max_dmov);
	wrapper_veo_args_set_float (kernel_ga_arg_ptr, narg++, dockpars.abs_max_dang);
	wrapper_veo_args_set_float (kernel_ga_arg_ptr, narg++, two_absmaxdmov);
	wrapper_veo_args_set_float (kernel_ga_arg_ptr, narg++, two_absmaxdang);
	wrapper_veo_args_set_float (kernel_ga_arg_ptr, narg++, dockpars.crossover_rate);
	wrapper_veo_args_set_u32   (kernel_ga_arg_ptr, narg++, dockpars.num_of_lsentities);
	wrapper_veo_args_set_u8    (kernel_ga_arg_ptr, narg++, dockpars.num_of_genes);
	
	// PC
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_pc_rotlist_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_pc_ref_coords_x_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_pc_ref_coords_y_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_pc_ref_coords_z_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_pc_rotbonds_moving_vectors_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_pc_rotbonds_unit_vectors_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_pc_ref_orientation_quats_const);
	wrapper_veo_args_set_u32 	(kernel_ga_arg_ptr, narg++, dockpars.rotbondlist_length);

	// IA
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_ie_atom_charges_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_ie_atom_types_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_contributors_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_reqm_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_reqm_hbond_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_atom1_types_reqm_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_atom2_types_reqm_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_VWpars_AC_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_VWpars_BD_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_dspars_S_const);
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_ia_dspars_V_const);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, dockpars.smooth);
	wrapper_veo_args_set_u32 	(kernel_ga_arg_ptr, narg++, dockpars.num_of_intraE_contributors);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, dockpars.grid_spacing);
	wrapper_veo_args_set_u32	(kernel_ga_arg_ptr, narg++, dockpars.num_of_atypes);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, dockpars.coeff_elec);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, dockpars.qasp);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, dockpars.coeff_desolv);

	// IE
	wrapper_veo_args_set_u64	(kernel_ga_arg_ptr, narg++, mem_dockpars_fgrids);
	wrapper_veo_args_set_u8     (kernel_ga_arg_ptr, narg++, dockpars.g1);
	wrapper_veo_args_set_u32 	(kernel_ga_arg_ptr, narg++, dockpars.g2);
	wrapper_veo_args_set_u32 	(kernel_ga_arg_ptr, narg++, dockpars.g3);
	wrapper_veo_args_set_u8     (kernel_ga_arg_ptr, narg++, dockpars.num_of_atoms);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, fgridsizex_minus1);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, fgridsizey_minus1);
	wrapper_veo_args_set_float 	(kernel_ga_arg_ptr, narg++, fgridsizez_minus1);
	wrapper_veo_args_set_u32 	(kernel_ga_arg_ptr, narg++, mul_tmp2);
	wrapper_veo_args_set_u32 	(kernel_ga_arg_ptr, narg++, mul_tmp3);

	// LS
	wrapper_veo_args_set_u16	(kernel_ga_arg_ptr, narg++, Host_max_num_of_iters);
	wrapper_veo_args_set_float	(kernel_ga_arg_ptr, narg++, dockpars.rho_lower_bound);
	wrapper_veo_args_set_float	(kernel_ga_arg_ptr, narg++, dockpars.base_dmov_mul_sqrt3);
	wrapper_veo_args_set_float	(kernel_ga_arg_ptr, narg++, dockpars.base_dang_mul_sqrt3);
	wrapper_veo_args_set_u8	    (kernel_ga_arg_ptr, narg++, Host_cons_limit);

	// Values changing every LGA run
	wrapper_veo_args_set_u32	(kernel_ga_arg_ptr, narg++, mypars->num_of_runs);


	// -----------------------------------------------------------------------------------------------------
	// Displaying kernel argument values
	// -----------------------------------------------------------------------------------------------------

	std::cout << "\n---------------------------------------------------------------------------------\n";
	std::cout << "Kernel LGA" << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n";
	std::cout << std::left << std::setw(40) << "mem_dockpars_conformations_current_Initial" << std::right << std::setw(31) << mem_dockpars_conformations_current_Initial << "\n";
	std::cout << std::left << std::setw(40) << "mem_dockpars_conformations_current_Final" << std::right << std::setw(31) << mem_dockpars_conformations_current_Final << "\n";
	std::cout << std::left << std::setw(40) << "mem_dockpars_energies_current" << std::right << std::setw(31) << mem_dockpars_energies_current << "\n";
	std::cout << std::left << std::setw(40) << "mem_evals_performed" << std::right << std::setw(31) << mem_evals_performed << "\n";
	std::cout << std::left << std::setw(40) << "mem_gens_performed" << std::right << std::setw(31) << mem_gens_performed << "\n";
	std::cout << std::left << std::setw(40) << "mem_prng_states" << std::right << std::setw(31) << mem_prng_states << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.pop_size" << std::right << std::setw(31) << dockpars.pop_size << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.num_of_energy_evals" << std::right << std::setw(31) << dockpars.num_of_energy_evals << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.num_of_generations" << std::right << std::setw(31) << dockpars.num_of_generations << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.tournament_rate" << std::right << std::setw(31) << dockpars.tournament_rate << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.mutation_rate" << std::right << std::setw(31) << dockpars.mutation_rate << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.abs_max_dmov" << std::right << std::setw(31) << dockpars.abs_max_dmov << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.abs_max_dang" << std::right << std::setw(31) << dockpars.abs_max_dang << "\n";
	std::cout << std::left << std::setw(40) << "two_absmaxdmov" << std::right << std::setw(31) << two_absmaxdmov << "\n";
	std::cout << std::left << std::setw(40) << "two_absmaxdang" << std::right << std::setw(31) << two_absmaxdang << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.crossover_rate" << std::right << std::setw(31) << dockpars.crossover_rate << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.num_of_lsentities" << std::right << std::setw(31) << dockpars.num_of_lsentities << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.num_of_genes" << std::right << std::setw(31) << int { dockpars.num_of_genes } << "\n";

	// PC
	std::cout << std::left << std::setw(40) << "mem_pc_rotlist_const" << std::right << std::setw(31) << mem_pc_rotlist_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_pc_ref_coords_x_const" << std::right << std::setw(31) << mem_pc_ref_coords_x_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_pc_ref_coords_y_const" << std::right << std::setw(31) << mem_pc_ref_coords_y_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_pc_ref_coords_z_const" << std::right << std::setw(31) << mem_pc_ref_coords_z_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_pc_rotbonds_moving_vectors_const" << std::right << std::setw(31) << mem_pc_rotbonds_moving_vectors_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_pc_rotbonds_unit_vectors_const" << std::right << std::setw(31) << mem_pc_rotbonds_unit_vectors_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_pc_ref_orientation_quats_const" << std::right << std::setw(31) << mem_pc_ref_orientation_quats_const << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.rotbondlist_length" << std::right << std::setw(31) << dockpars.rotbondlist_length << "\n";

	// IA
	std::cout << std::left << std::setw(40) << "mem_ia_ie_atom_charges_const" << std::right << std::setw(31) << mem_ia_ie_atom_charges_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_ie_atom_types_const" << std::right << std::setw(31) << mem_ia_ie_atom_types_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_contributors_const" << std::right << std::setw(31) << mem_ia_contributors_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_reqm_const" << std::right << std::setw(31) << mem_ia_reqm_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_reqm_hbond_const" << std::right << std::setw(31) << mem_ia_reqm_hbond_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_atom1_types_reqm_const" << std::right << std::setw(31) << mem_ia_atom1_types_reqm_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_atom2_types_reqm_const" << std::right << std::setw(31) << mem_ia_atom2_types_reqm_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_VWpars_AC_const" << std::right << std::setw(31) << mem_ia_VWpars_AC_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_VWpars_BD_const" << std::right << std::setw(31) << mem_ia_VWpars_BD_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_dspars_S_const" << std::right << std::setw(31) << mem_ia_dspars_S_const << "\n";
	std::cout << std::left << std::setw(40) << "mem_ia_dspars_V_const" << std::right << std::setw(31) << mem_ia_dspars_V_const << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.smooth" << std::right << std::setw(31) << dockpars.smooth << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.num_of_intraE_contributors" << std::right << std::setw(31) << dockpars.num_of_intraE_contributors << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.grid_spacing" << std::right << std::setw(31) << dockpars.grid_spacing << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.num_of_atypes" << std::right << std::setw(31) << dockpars.num_of_atypes << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.coeff_elec" << std::right << std::setw(31) << dockpars.coeff_elec << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.qasp" << std::right << std::setw(31) << dockpars.qasp << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.coeff_desolv" << std::right << std::setw(31) << dockpars.coeff_desolv << "\n";

	// IE
	std::cout << std::left << std::setw(40) << "mem_dockpars_fgrids" << std::right << std::setw(31) << mem_dockpars_fgrids << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.g1" << std::right << std::setw(31) << int { dockpars.g1 } << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.g2" << std::right << std::setw(31) << dockpars.g2 << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.g3" << std::right << std::setw(31) << dockpars.g3 << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.num_of_atoms" << std::right << std::setw(31) << int { dockpars.num_of_atoms } << "\n";
	std::cout << std::left << std::setw(40) << "fgridsizex_minus1" << std::right << std::setw(31) << fgridsizex_minus1 << "\n";
	std::cout << std::left << std::setw(40) << "fgridsizey_minus1" << std::right << std::setw(31) << fgridsizey_minus1 << "\n";
	std::cout << std::left << std::setw(40) << "fgridsizez_minus1" << std::right << std::setw(31) << fgridsizez_minus1 << "\n";
	std::cout << std::left << std::setw(40) << "mul_tmp2" << std::right << std::setw(31) << mul_tmp2 << "\n";
	std::cout << std::left << std::setw(40) << "mul_tmp3" << std::right << std::setw(31) << mul_tmp3 << "\n";

	// LS
	std::cout << std::left << std::setw(40) << "Host_max_num_of_iters" << std::right << std::setw(31) << Host_max_num_of_iters << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.rho_lower_bound" << std::right << std::setw(31) << dockpars.rho_lower_bound << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.base_dmov_mul_sqrt3" << std::right << std::setw(31) << dockpars.base_dmov_mul_sqrt3 << "\n";
	std::cout << std::left << std::setw(40) << "dockpars.base_dang_mul_sqrt3" << std::right << std::setw(31) << dockpars.base_dang_mul_sqrt3 << "\n";
	std::cout << std::left << std::setw(40) << "Host_cons_limit" << std::right << std::setw(31) << int { Host_cons_limit } << "\n";
	std::cout << std::left << std::setw(40) << "mypars->num_of_runs" << std::right << std::setw(31) << mypars->num_of_runs << "\n";
	std::cout << std::endl;
	std::cout << "\n---------------------------------------------------------------------------------\n";

	// -----------------------------------------------------------------------------------------------------
	// Launching kernel
	// -----------------------------------------------------------------------------------------------------

	uint64_t kernel_ga_id;
	uint64_t retval_ga;

	clock_start_docking = clock();

	printf("Docking runs to be executed: %lu\n", mypars->num_of_runs); 
	printf("Execution run: ");

	kernel_ga_id = wrapper_veo_call_async_by_name(veo_thread_context, kernel_ga_handle, name_k_ga, kernel_ga_arg_ptr);
	wrapper_veo_call_wait_result(veo_thread_context, kernel_ga_id, &retval_ga);

	clock_stop_docking = clock();

	printf("\n");
	fflush(stdout);

	// -----------------------------------------------------------------------------------------------------
	// Reading results from device
	// -----------------------------------------------------------------------------------------------------

	wrapper_veo_read_mem (ve_process, cpu_final_populations.data(), mem_dockpars_conformations_current_Final, size_populations_nbytes);
	wrapper_veo_read_mem (ve_process, cpu_energies.data(), mem_dockpars_energies_current, size_energies_nbytes);
	wrapper_veo_read_mem (ve_process, cpu_evals_of_runs.data(), mem_evals_performed, size_evals_of_runs_nbytes);
	wrapper_veo_read_mem (ve_process, cpu_gens_of_runs.data(), mem_gens_performed, size_evals_of_runs_nbytes);

	// -----------------------------------------------------------------------------------------------------
	// Processing results
	// -----------------------------------------------------------------------------------------------------

	/*
	for (int cnt_pop = 0; cnt_pop < size_populations_nelems; cnt_pop++)
		printf("total_num_pop: %u, cpu_final_populations[%u]: %f\n", size_populations_nelems, cnt_pop, cpu_final_populations[cnt_pop]);

	for (int cnt_pop = 0; cnt_pop < size_energies_nelems; cnt_pop++)
		printf("total_num_energies: %u, cpu_energies[%u]: %f\n", size_energies_nelems, cnt_pop, cpu_energies[cnt_pop]);
	*/

	for (unsigned int run_cnt = 0; run_cnt < mypars->num_of_runs; run_cnt++) {

		arrange_result(
			cpu_final_populations.data() + run_cnt * mypars->pop_size * ACTUAL_GENOTYPE_LENGTH,
			cpu_energies.data()          + run_cnt * mypars->pop_size,
			mypars->pop_size);

		/*printf("cpu_evals_of_runs[%u]: %u\n", run_cnt, cpu_evals_of_runs[run_cnt]);*/

		make_resfiles(
			cpu_final_populations.data() + run_cnt * mypars->pop_size * ACTUAL_GENOTYPE_LENGTH,
			cpu_energies.data()          + run_cnt * mypars->pop_size,
			&myligand_reference,
			myligand_init,
			myxrayligand,
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
	for (int cnt_pop = 0; cnt_pop < size_populations_nelems; cnt_pop++)
		printf("total_num_pop: %u, cpu_final_populations[%u]: %f\n", size_populations_nelems, cnt_pop, cpu_final_populations[cnt_pop]);

	for (int cnt_pop = 0; cnt_pop < size_energies_nelems; cnt_pop++)
		printf("total_num_energies: %u, cpu_energies[%u]: %f\n", size_energies_nelems, cnt_pop, cpu_energies[cnt_pop]);
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
