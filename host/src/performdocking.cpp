#include "performdocking.h"
#include "device_args.h"
#include "packbuff.hpp"
#include <stddef.h>
#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#define STRINGIZE2(s) #s
#define STRINGIZE(s)	STRINGIZE2(s)
#define KRNL_BIN_FOLDER STRINGIZE(KRNL_LIB_DIRECTORY)
#define KRNL_SRC_FOLDER STRINGIZE(KRNL_SRC_DIRECTORY)
#define KRNL_SRC_COMMON STRINGIZE(KCMN_SRC_DIRECTORY)
#define KRNL0 	STRINGIZE(K0)
#define KRNL_GA STRINGIZE(K_GA)

// Spaces for std::cout
#define SPACE_L	50
#define SPACE_M	31
#define SPACE_S	10

// Docking
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
	// -------------------------------------------------------------------------
	// Host Setup
	// -------------------------------------------------------------------------

	const char* name_k_ga = KRNL_GA;
	const char* krnl_folder = KRNL_BIN_FOLDER;
	char path_k_ga[100];
	strcpy(path_k_ga, krnl_folder);
	strcat(path_k_ga, "/");
	strcat(path_k_ga, name_k_ga);
	strcat(path_k_ga, ".so");
	std::cout << std::endl;
	std::cout << std::endl << "Kernel: " << path_k_ga << std::endl;
	std::cout << std::endl;

	// Determine number of VE processes to start
	// If starting multiple on same VE, OpenMP thread number must be set appropriately
	vector<int> ve_node_ids;

	if (const char* env_p = std::getenv("VE_NODE_IDS")) {
		std::stringstream env_s;
		env_s << std::string(env_p);
		while (env_s.good()) {
			std::string substr;
			getline(env_s, substr, ',');
			ve_node_ids.push_back(std::stoi(substr));
		}
	} else {
		ve_node_ids.push_back(0);
	}
	std::cout << "Running in " << ve_node_ids.size() << " VE processes" << std::endl;
	std::cout << "VE numbers ";
	for (int v : ve_node_ids)
		std::cout << v << " ";
	std::cout << std::endl;

	int ve_num_procs = ve_node_ids.size();

	// VEO code
	wrapper_veo_api_version ();

	// Start VE procs
	veo_proc_handle *ve_process[ve_num_procs];
	uint64_t kernel_ga_handle[ve_num_procs];
	veo_thr_ctxt *veo_thread_context[ve_num_procs];

	for (int i = 0; i < ve_num_procs; i++) {
		ve_process[i] = wrapper_veo_proc_create(ve_node_ids[i]);
		kernel_ga_handle[i] = wrapper_veo_load_library(ve_process[i], path_k_ga);
		veo_thread_context[i] = wrapper_veo_context_open(ve_process[i]);
	}

	// Maximum number of runs per VE proc
	int ve_num_of_runs = (mypars->num_of_runs + ve_num_procs - 1) / ve_num_procs;

	// End of Host Setup
	// -------------------------------------------------------------------------

	clock_t clock_start_docking;
	clock_t	clock_stop_docking;
	clock_t clock_stop_program_before_clustering;

	timeval tv_start_docking, tv_stop_docking;

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
	size_t size_prng_seeds_nelems = mypars->pop_size * ve_num_procs;
	size_t size_prng_seeds_nbytes = size_prng_seeds_nelems * sizeof(unsigned int);
	std::vector<unsigned int> cpu_prng_seeds (size_prng_seeds_nelems);
	
	// Initializing seed generator
	genseed(time(NULL));

	// Generating seeds (for each thread during GA)
	for (int ve_id = 0; ve_id < ve_num_procs; ve_id++) {
#ifdef REPRO
		// each VE process gets the same seeds in reproducibility mode
		// this only makes sense for testing and debugging!
		genseed(1234567u);
#endif
		for (int i = 0; i < mypars->pop_size; i++) {
			cpu_prng_seeds[ve_id * mypars->pop_size + i] = genseed(0u);
		}
	}

	size_t size_evals_of_runs_nelems = mypars->num_of_runs;
	size_t size_evals_of_runs_nbytes = size_evals_of_runs_nelems * sizeof(int);

	// Allocating memory in CPU for evaluation counters
	std::vector<int> cpu_evals_of_runs (size_evals_of_runs_nelems, 0);
	
	// Allocating memory in CPU for generation counters
	std::vector<int> cpu_gens_of_runs (size_evals_of_runs_nelems, 0);

	// -------------------------------------------------------------------------
	// Preparing the constant data fields for the accelerator (calcenergy.cpp)
	// -------------------------------------------------------------------------

	kernelconstant_static  KerConstStatic;
	if (prepare_conststatic_fields_for_aurora(&myligand_reference, mypars, cpu_ref_ori_angles.data(), &KerConstStatic) == 1)
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

	// -------------------------------------------------------------------------
	// Defining kernel buffers
	// -------------------------------------------------------------------------

	// LGA

	// Pose Calculation buffers
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

	// IE buffers
	size_t size_InterE_atom_charges_nelems = MAX_NUM_OF_ATOMS;
	size_t size_InterE_atom_charges_nbytes = size_InterE_atom_charges_nelems * sizeof(float);

	size_t size_InterE_atom_types_nelems = MAX_NUM_OF_ATOMS;
	size_t size_InterE_atom_types_nbytes = size_InterE_atom_types_nelems * sizeof(char);

	// IA buffers
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

	// -------------------------------------------------------------------------
	// Printing sizes
	// -------------------------------------------------------------------------

	// LGA
#ifdef DOCK_DEBUG
	std::cout << "\n---------------------------------------------------------------------------------\n";
	std::cout << std::left << std::setw(SPACE_L) << "Memory sizes" << std::right << std::setw(SPACE_M) << "Bytes" << std::right << std::setw(SPACE_S) << "KB" << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n";
	std::cout << std::left << std::setw(SPACE_L) << "size_populations_nbytes" << std::right << std::setw(SPACE_M) << size_populations_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_populations_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_energies_nbytes" << std::right << std::setw(SPACE_M) << size_energies_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_energies_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_evals_of_runs_nbytes" << std::right << std::setw(SPACE_M) << size_evals_of_runs_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_evals_of_runs_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_prng_seeds_nbytes" << std::right << std::setw(SPACE_M) << size_prng_seeds_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_prng_seeds_nbytes) << std::endl;

	size_t size_ga = size_populations_nbytes + size_populations_nbytes +
		size_energies_nbytes + size_evals_of_runs_nbytes +
		size_evals_of_runs_nbytes + size_prng_seeds_nbytes;

	// Pose Calculation
	std::cout << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_rotlist_nbytes" << std::right << std::setw(SPACE_M) << size_rotlist_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_rotlist_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_ref_coords_nbytes" << std::right << std::setw(SPACE_M) << size_ref_coords_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_ref_coords_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_rotbonds_moving_vectors_nbytes" << std::right << std::setw(SPACE_M) << size_rotbonds_moving_vectors_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_rotbonds_moving_vectors_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_rotbonds_unit_vectors_nbytes" << std::right << std::setw(SPACE_M) << size_rotbonds_unit_vectors_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_rotbonds_unit_vectors_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_ref_orientation_quats_nbytes" << std::right << std::setw(SPACE_M) << size_ref_orientation_quats_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_ref_orientation_quats_nbytes) << std::endl;

	size_t size_pc = size_rotlist_nbytes + size_ref_coords_nbytes +
		size_ref_coords_nbytes + size_ref_coords_nbytes +
		size_rotbonds_moving_vectors_nbytes + size_rotbonds_unit_vectors_nbytes +
		size_ref_orientation_quats_nbytes;

	// IA
	std::cout << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_InterE_atom_charges_nbytes" << std::right << std::setw(SPACE_M) << size_InterE_atom_charges_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_InterE_atom_charges_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_InterE_atom_types_nbytes" << std::right << std::setw(SPACE_M) << size_InterE_atom_types_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_InterE_atom_types_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_intraE_contributors_nbytes" << std::right << std::setw(SPACE_M) << size_intraE_contributors_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_intraE_contributors_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_reqm_nbytes" << std::right << std::setw(SPACE_M) << size_reqm_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_reqm_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_reqm_hbond_nbytes" << std::right << std::setw(SPACE_M) << size_reqm_hbond_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_reqm_hbond_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_atom1_types_reqm_nbytes" << std::right << std::setw(SPACE_M) << size_atom1_types_reqm_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_atom1_types_reqm_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_atom2_types_reqm_nbytes" << std::right << std::setw(SPACE_M) << size_atom2_types_reqm_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_atom2_types_reqm_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_VWpars_AC_nbytes" << std::right << std::setw(SPACE_M) << size_VWpars_AC_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_VWpars_AC_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_VWpars_BD_nbytes" << std::right << std::setw(SPACE_M) << size_VWpars_BD_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_VWpars_BD_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_dspars_S_nbytes" << std::right << std::setw(SPACE_M) << size_dspars_S_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_dspars_S_nbytes) << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_dspars_V_nbytes" << std::right << std::setw(SPACE_M) << size_dspars_V_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_dspars_V_nbytes) << std::endl;

	size_t size_ia = size_InterE_atom_charges_nbytes + size_InterE_atom_types_nbytes +
		size_intraE_contributors_nbytes + size_reqm_nbytes + size_reqm_hbond_nbytes +
		size_atom1_types_reqm_nbytes + size_atom2_types_reqm_nbytes +
		size_VWpars_AC_nbytes + size_VWpars_BD_nbytes +
		size_dspars_S_nbytes + size_dspars_V_nbytes;

	// IE
	std::cout << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "size_floatgrids_nbytes" << std::right << std::setw(SPACE_M) << size_floatgrids_nbytes << std::right << std::setw(SPACE_S) << sizeKB(size_floatgrids_nbytes) << std::endl;

	size_t size_ie = size_floatgrids_nbytes;

	// Total amount memory
	size_t size_total_mem = size_ga + size_pc + size_ia + size_ie;

	std::cout << std::endl;
	std::cout << std::left << std::setw(SPACE_L) << "Total memory" << std::right << std::setw(SPACE_M) << size_total_mem << std::right << std::setw(SPACE_S) << sizeKB(size_total_mem) << std::endl;
	std::cout << "---------------------------------------------------------------------------------\n" << std::endl;
#endif

	// -------------------------------------------------------------------------
	// Defining kernel arguments
	// -------------------------------------------------------------------------

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

	// The device_args structure contains the arguments for the kernel function
	device_args da;
	// VE virtual address of da structure, passed to kernel
	uint64_t da_VEMVA;

	// VE proc ID, is changed for each new proc in multi-proc runs
	da.ve_proc_id = 0;
	da.ve_num_procs = ve_num_procs;

	// Filling scalar values into device_args structure
	// GA
	da.DockConst_pop_size            = dockpars.pop_size;
	da.DockConst_num_of_energy_evals = dockpars.num_of_energy_evals;
	da.DockConst_num_of_generations  = dockpars.num_of_generations;
	da.DockConst_tournament_rate     = dockpars.tournament_rate;
	da.DockConst_mutation_rate       = dockpars.mutation_rate;
	da.DockConst_abs_max_dmov        = dockpars.abs_max_dmov;    
	da.DockConst_abs_max_dang        = dockpars.abs_max_dang;
	da.Host_two_absmaxdmov           = two_absmaxdmov;
	da.Host_two_absmaxdang           = two_absmaxdang;
	da.DockConst_crossover_rate      = dockpars.crossover_rate;
	da.DockConst_num_of_lsentities   = dockpars.num_of_lsentities;
	da.DockConst_num_of_genes        = dockpars.num_of_genes;
	// PC
	da.DockConst_rotbondlist_length  = dockpars.rotbondlist_length;
	// IA
	da.DockConst_smooth                     = dockpars.smooth;
	da.DockConst_num_of_intraE_contributors = dockpars.num_of_intraE_contributors;
	da.DockConst_grid_spacing               = dockpars.grid_spacing;
	da.DockConst_num_of_atypes              = dockpars.num_of_atypes;
	da.DockConst_coeff_elec                 = dockpars.coeff_elec;
	da.DockConst_qasp                       = dockpars.qasp;
	da.DockConst_coeff_desolv               = dockpars.coeff_desolv;
	// IE
	da.DockConst_xsz               = dockpars.gridsize_x;
	da.DockConst_ysz               = dockpars.gridsize_y;
	da.DockConst_zsz               = dockpars.gridsize_z;
	da.DockConst_num_of_atoms      = dockpars.num_of_atoms;
	da.DockConst_gridsize_x_minus1 = fgridsizex_minus1;
	da.DockConst_gridsize_y_minus1 = fgridsizey_minus1;
	da.DockConst_gridsize_z_minus1 = fgridsizez_minus1;
	da.Host_mul_tmp2               = mul_tmp2;
	da.Host_mul_tmp3               = mul_tmp3;
	// LS
	da.DockConst_max_num_of_iters    = Host_max_num_of_iters;
	da.DockConst_rho_lower_bound     = dockpars.rho_lower_bound; 
	da.DockConst_base_dmov_mul_sqrt3 = dockpars.base_dmov_mul_sqrt3;
	da.DockConst_base_dang_mul_sqrt3 = dockpars.base_dang_mul_sqrt3;
	da.DockConst_cons_limit          = Host_cons_limit;
	// Values changing every LGA run
	da.Host_num_of_runs = mypars->num_of_runs;
	
	// Packed buffer with all arguments to be passed, alloc chunksize = 32MB
	auto pb = PackBuff(32*1024*1024);
	
	// Pack the device_args struct as first element into the packed buffer.
	// This way we know how to find it easily
	pb.pack((void *)&da, sizeof(da), (uint64_t)&da_VEMVA);

#define DA_IN_PB_PTR(ELEMENT) ((uint64_t)pb.data() + offsetof(device_args, ELEMENT))

	// GA
	pb.pack(cpu_init_populations.data(), size_populations_nbytes, DA_IN_PB_PTR(PopulationCurrentInitial));
	pb.pack(NULL, size_populations_nbytes, DA_IN_PB_PTR(PopulationCurrentFinal));
	pb.pack(NULL, size_energies_nbytes, DA_IN_PB_PTR(EnergyCurrent));
	pb.pack(NULL, size_evals_of_runs_nbytes, DA_IN_PB_PTR(Evals_performed));
	pb.pack(NULL, size_evals_of_runs_nbytes, DA_IN_PB_PTR(Gens_performed));
	pb.pack(cpu_prng_seeds.data(), size_prng_seeds_nbytes, DA_IN_PB_PTR(dockpars_prng_states));
	
	// PC
	pb.pack(&KerConstStatic.rotlist_const[0], size_rotlist_nbytes, DA_IN_PB_PTR(PC_rotlist));
	pb.pack(&KerConstStatic.ref_coords_x_const[0], size_ref_coords_nbytes, DA_IN_PB_PTR(PC_ref_coords_x));
	pb.pack(&KerConstStatic.ref_coords_y_const[0], size_ref_coords_nbytes, DA_IN_PB_PTR(PC_ref_coords_y));
	pb.pack(&KerConstStatic.ref_coords_z_const[0], size_ref_coords_nbytes, DA_IN_PB_PTR(PC_ref_coords_z));
	pb.pack(&KerConstStatic.rotbonds_moving_vectors_const[0], size_rotbonds_moving_vectors_nbytes, DA_IN_PB_PTR(PC_rotbonds_moving_vectors));
	pb.pack(&KerConstStatic.rotbonds_unit_vectors_const[0], size_rotbonds_unit_vectors_nbytes, DA_IN_PB_PTR(PC_rotbonds_unit_vectors));
	pb.pack(&KerConstStatic.ref_orientation_quats_const[0], size_ref_orientation_quats_nbytes, DA_IN_PB_PTR(PC_ref_orientation_quats));

	// IA
	pb.pack(&KerConstStatic.atom_charges_const[0], size_InterE_atom_charges_nbytes, DA_IN_PB_PTR(IA_IE_atom_charges));
	pb.pack(&KerConstStatic.atom_types_const[0], size_InterE_atom_types_nbytes, DA_IN_PB_PTR(IA_IE_atom_types));
	pb.pack(&KerConstStatic.intraE_contributors_const[0], size_intraE_contributors_nbytes, DA_IN_PB_PTR(IA_intraE_contributors));
	pb.pack(&KerConstStatic.reqm_const[0], size_reqm_nbytes, DA_IN_PB_PTR(IA_reqm));
	pb.pack(&KerConstStatic.reqm_hbond_const[0], size_reqm_hbond_nbytes, DA_IN_PB_PTR(IA_reqm_hbond));
	pb.pack(&KerConstStatic.atom1_types_reqm_const[0], size_atom1_types_reqm_nbytes, DA_IN_PB_PTR(IA_atom1_types_reqm));
	pb.pack(&KerConstStatic.atom2_types_reqm_const[0], size_atom2_types_reqm_nbytes, DA_IN_PB_PTR(IA_atom2_types_reqm));
	pb.pack(&KerConstStatic.VWpars_AC_const[0], size_VWpars_AC_nbytes, DA_IN_PB_PTR(IA_VWpars_AC));
	pb.pack(&KerConstStatic.VWpars_BD_const[0], size_VWpars_BD_nbytes, DA_IN_PB_PTR(IA_VWpars_BD));
	pb.pack(&KerConstStatic.dspars_S_const[0], size_dspars_S_nbytes, DA_IN_PB_PTR(IA_dspars_S));
	pb.pack(&KerConstStatic.dspars_V_const[0], size_dspars_V_nbytes, DA_IN_PB_PTR(IA_dspars_V));

	// IE
	pb.pack(cpu_floatgrids, size_floatgrids_nbytes, DA_IN_PB_PTR(Fgrids));

	// We keep a copy of each proc's device_args structure
	device_args *da_copy[ve_num_procs];
	
	uint64_t kernel_ga_id[ve_num_procs];
	uint64_t retval_ga[ve_num_procs];
	uint64_t padding[ve_num_procs];
#ifndef PADDING
#define PADDING 0
#endif
	for (int ve_id = 0; ve_id < ve_num_procs; ve_id++) {
		// Allocating memory for entire packed buffer on VE
		uint64_t mem_kernel_args_packbuff;
		wrapper_veo_alloc_mem(ve_process[ve_id], &mem_kernel_args_packbuff, pb.size() + 128 * ve_num_procs);

		//// save data buffer
		//pb.save("packbuff_save.dat");

		// Setting VE proc ID in packbuff
		device_args *pb_da = (device_args *)pb.data();
		pb_da->ve_proc_id = ve_id;

		// Padding
		padding[ve_id] = PADDING * ve_id;
		
		// Fixing "relocation" addresses to point to VE virtual addresses
		pb.fixup(mem_kernel_args_packbuff + padding[ve_id]);

		// Updating local device_args structure (copy) such that we can free the packbuff later
		da_copy[ve_id] = new device_args;
		memcpy(da_copy[ve_id], pb.data(), sizeof(da));

		// Transferring packbuff to device
		wrapper_veo_write_mem(ve_process[ve_id], mem_kernel_args_packbuff + padding[ve_id], pb.data(), pb.size());

		// TODO: this is a memory leak because there is no free
		veo_args *kernel_ga_arg_ptr = wrapper_veo_args_alloc();

		// Creating a VEO arguments object
		wrapper_veo_args_set_u64(kernel_ga_arg_ptr, 0, mem_kernel_args_packbuff + padding[ve_id]);

#ifdef DOCK_DEBUG
		if (ve_id == 0) {
			// -----------------------------------------------------------------
			// Displaying kernel argument values
			// -----------------------------------------------------------------

			std::cout << "\n---------------------------------------------------------------------------------\n";
			std::cout << "Kernel LGA" << std::endl;
			std::cout << "---------------------------------------------------------------------------------\n";

			std::cout << std::left << std::setw(SPACE_L) << "mem_dockpars_conformations_current_Final" << std::right << std::setw(SPACE_M) << da.PopulationCurrentFinal << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "mem_dockpars_energies_current" << std::right << std::setw(SPACE_M) << da.EnergyCurrent << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "mem_evals_performed" << std::right << std::setw(SPACE_M) << da.Evals_performed << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "mem_gens_performed" << std::right << std::setw(SPACE_M) << da.Gens_performed << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.pop_size" << std::right << std::setw(SPACE_M) << dockpars.pop_size << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.num_of_energy_evals" << std::right << std::setw(SPACE_M) << dockpars.num_of_energy_evals << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.num_of_generations" << std::right << std::setw(SPACE_M) << dockpars.num_of_generations << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.tournament_rate" << std::right << std::setw(SPACE_M) << dockpars.tournament_rate << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.mutation_rate" << std::right << std::setw(SPACE_M) << dockpars.mutation_rate << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.abs_max_dmov" << std::right << std::setw(SPACE_M) << dockpars.abs_max_dmov << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.abs_max_dang" << std::right << std::setw(SPACE_M) << dockpars.abs_max_dang << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "two_absmaxdmov" << std::right << std::setw(SPACE_M) << two_absmaxdmov << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "two_absmaxdang" << std::right << std::setw(SPACE_M) << two_absmaxdang << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.crossover_rate" << std::right << std::setw(SPACE_M) << dockpars.crossover_rate << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.num_of_lsentities" << std::right << std::setw(SPACE_M) << dockpars.num_of_lsentities << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.num_of_genes" << std::right << std::setw(SPACE_M) << int { dockpars.num_of_genes } << "\n";
			std::cout << std::endl;

			// PC
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.rotbondlist_length" << std::right << std::setw(SPACE_M) << dockpars.rotbondlist_length << "\n";
			std::cout << std::endl;

			// IA
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.smooth" << std::right << std::setw(SPACE_M) << dockpars.smooth << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.num_of_intraE_contributors" << std::right << std::setw(SPACE_M) << dockpars.num_of_intraE_contributors << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.grid_spacing" << std::right << std::setw(SPACE_M) << dockpars.grid_spacing << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.num_of_atypes" << std::right << std::setw(SPACE_M) << dockpars.num_of_atypes << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.coeff_elec" << std::right << std::setw(SPACE_M) << dockpars.coeff_elec << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.qasp" << std::right << std::setw(SPACE_M) << dockpars.qasp << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.coeff_desolv" << std::right << std::setw(SPACE_M) << dockpars.coeff_desolv << "\n";
			std::cout << std::endl;

			// IE
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.g1" << std::right << std::setw(SPACE_M) << int { dockpars.g1 } << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.g2" << std::right << std::setw(SPACE_M) << dockpars.g2 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.g3" << std::right << std::setw(SPACE_M) << dockpars.g3 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.num_of_atoms" << std::right << std::setw(SPACE_M) << int { dockpars.num_of_atoms } << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "fgridsizex_minus1" << std::right << std::setw(SPACE_M) << fgridsizex_minus1 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "fgridsizey_minus1" << std::right << std::setw(SPACE_M) << fgridsizey_minus1 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "fgridsizez_minus1" << std::right << std::setw(SPACE_M) << fgridsizez_minus1 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "mul_tmp2" << std::right << std::setw(SPACE_M) << mul_tmp2 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "mul_tmp3" << std::right << std::setw(SPACE_M) << mul_tmp3 << "\n";
			std::cout << std::endl;

			// LS
			std::cout << std::left << std::setw(SPACE_L) << "Host_max_num_of_iters" << std::right << std::setw(SPACE_M) << Host_max_num_of_iters << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.rho_lower_bound" << std::right << std::setw(SPACE_M) << dockpars.rho_lower_bound << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.base_dmov_mul_sqrt3" << std::right << std::setw(SPACE_M) << dockpars.base_dmov_mul_sqrt3 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "dockpars.base_dang_mul_sqrt3" << std::right << std::setw(SPACE_M) << dockpars.base_dang_mul_sqrt3 << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "Host_cons_limit" << std::right << std::setw(SPACE_M) << int { Host_cons_limit } << "\n";
			std::cout << std::left << std::setw(SPACE_L) << "mypars->num_of_runs" << std::right << std::setw(SPACE_M) << mypars->num_of_runs << "\n";
			std::cout << "\n---------------------------------------------------------------------------------\n";
			std::cout << std::endl;
		}
#endif

		// ---------------------------------------------------------------------
		// Launching kernel
		// ---------------------------------------------------------------------

		if (ve_id == 0) {
			clock_start_docking = clock();
			gettimeofday(&tv_start_docking, NULL);
		}

		std::cout << std::endl;
		std::cout << "Docking runs to be executed: " << mypars->num_of_runs << std::endl;
		std::cout << "Execution run: (VE proc " << ve_id << " )" << std::endl;
		std::cout << std::flush;

		kernel_ga_id[ve_id] = wrapper_veo_call_async_by_name(veo_thread_context[ve_id], kernel_ga_handle[ve_id], name_k_ga, kernel_ga_arg_ptr);

		// SERIALIZE! for testing, only
		//wrapper_veo_call_wait_result(veo_thread_context[ve_id], kernel_ga_id[ve_id], &retval_ga[ve_id]);
	}
#if 1
	for (int ve_id = 0; ve_id < ve_num_procs; ve_id++) {
		wrapper_veo_call_wait_result(veo_thread_context[ve_id], kernel_ga_id[ve_id], &retval_ga[ve_id]);
	}
#endif

	clock_stop_docking = clock();
	gettimeofday(&tv_stop_docking, NULL);
	double elapsed_time_docking;
	elapsed_time_docking = (((double)tv_stop_docking.tv_sec * 1.e6 + (double)tv_stop_docking.tv_usec) -
				((double)tv_start_docking.tv_sec * 1.e6 + (double)tv_start_docking.tv_usec)) / 1.e6;

	std::cout << std::endl;
	std::cout << std::flush;

	// -------------------------------------------------------------------------
	// Reading results from device
	// -------------------------------------------------------------------------

	//////////
	// TODO EF
	//////////
	for (int ve_id = 0; ve_id < ve_num_procs; ve_id++) {
		int run_begin = ve_id * ve_num_of_runs;
		int run_end = std::min((ve_id + 1) * ve_num_of_runs, (int)mypars->num_of_runs);
		
		size_t size_populations_nelems = (run_end - run_begin) * mypars->pop_size * ACTUAL_GENOTYPE_LENGTH;
		size_t size_populations_nbytes = size_populations_nelems * sizeof(float);
		size_t size_populations_offs = run_begin * mypars->pop_size * ACTUAL_GENOTYPE_LENGTH * sizeof(float);
		size_t size_energies_nelems = (run_end - run_begin) * mypars->pop_size;
		size_t size_energies_nbytes = size_energies_nelems * sizeof(float);
		size_t size_energies_offs = run_begin * mypars->pop_size * sizeof(float);
		size_t size_evals_of_runs_nelems = run_end - run_begin;
		size_t size_evals_of_runs_nbytes = size_evals_of_runs_nelems * sizeof(int);
		size_t size_evals_of_runs_offs = run_begin * sizeof(int);

		wrapper_veo_read_mem (ve_process[ve_id], (void *)((char *)cpu_final_populations.data() + size_populations_offs),
				      (uint64_t)da_copy[ve_id]->PopulationCurrentFinal + size_populations_offs, size_populations_nbytes);
		wrapper_veo_read_mem (ve_process[ve_id], (void *)((char *)cpu_energies.data() + size_energies_offs),
				      (uint64_t)da_copy[ve_id]->EnergyCurrent + size_energies_offs, size_energies_nbytes);
		wrapper_veo_read_mem (ve_process[ve_id], (void *)((char *)cpu_evals_of_runs.data() + size_evals_of_runs_offs),
				      (uint64_t)da_copy[ve_id]->Evals_performed + size_evals_of_runs_offs, size_evals_of_runs_nbytes);
		wrapper_veo_read_mem (ve_process[ve_id], (void *)((char *)cpu_gens_of_runs.data() + size_evals_of_runs_offs),
				      (uint64_t)da_copy[ve_id]->Gens_performed + size_evals_of_runs_offs, size_evals_of_runs_nbytes);
	}

	// Destroying the VEO process early in order to get a more accurate PROGINF
	for (int ve_id = 0; ve_id < ve_num_procs; ve_id++) {
		veo_proc_destroy(ve_process[ve_id]);
	}

	// -------------------------------------------------------------------------
	// Processing results
	// -------------------------------------------------------------------------

	/*
	for (int cnt_pop = 0; cnt_pop < size_populations_nelems; cnt_pop++)
		printf("total_num_pop: %u, cpu_final_populations[%u]: %f\n", size_populations_nelems, cnt_pop, cpu_final_populations[cnt_pop]);

	for (int cnt_pop = 0; cnt_pop < size_energies_nelems; cnt_pop++)
		printf("total_num_energies: %u, cpu_energies[%u]: %f\n", size_energies_nelems, cnt_pop, cpu_energies[cnt_pop]);
	*/
	long sum_energy_evals = 0, sum_generations = 0;
	for (unsigned int run_cnt = 0; run_cnt < mypars->num_of_runs; run_cnt++) {
		sum_energy_evals += cpu_evals_of_runs[run_cnt];
		sum_generations += cpu_gens_of_runs[run_cnt];
	}
	printf("Time spent in docking search      : %.3fs\n", elapsed_time_docking);
	printf("Total number of energy evaluations: %ld -> %.5f us/eval\n", sum_energy_evals, (1.e6 * elapsed_time_docking) / (float)sum_energy_evals);
	printf("Total number of generations       : %ld\n", sum_generations);

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

        printf("arrange_result and make_resfiles took : %.3fs\n",
               ELAPSEDSECS(clock_stop_program_before_clustering, clock_stop_docking));

	clusanal_gendlg(cpu_result_ligands.data(), 
			mypars->num_of_runs,
			myligand_init, mypars,
   		        mygrid, 
			argc,
			argv,
			elapsed_time_docking,
			ELAPSEDSECS(clock_stop_program_before_clustering, clock_start_program));

	clock_stop_docking = clock();

        printf("clustering took : %.3fs\n",
               ELAPSEDSECS(clock_stop_docking, clock_stop_program_before_clustering));

	return 0;
}
