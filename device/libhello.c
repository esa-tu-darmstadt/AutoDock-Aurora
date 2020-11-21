#include <stdio.h>
#include <stdint.h>
#include <ftrace.h>

#include "defines.h"
/*
#include "libkernel1.c"
#include "libkernel2.c"
#include "libkernel4.c"
#include "libkernel3.c"
*/

uint64_t hello()
{
  printf("Hello, world\n");
  return 0;
}

int main() {
  printf("This is main()\n");
  hello();

/*
  // Using real VE VMA addresses for artifially invoking kernel2
  
  // Defining and initializing structs
  // as local variables -> will be stored in the stack!
  // This mimic what we did in our VEOffload design
  // as we passed some into the stack.
  // Initilization eference: host/src/calcenergy.cpp
  int m;

  kernelconstant_interintra K_interintra;
  for (m = 0; m < MAX_NUM_OF_ATOMS; m++) { 
    K_interintra.atom_charges_const[m] = (float) 0.1f; 
    K_interintra.atom_types_const[m] = (char) 5;
  }

  kernelconstant_intracontrib K_intracontrib;
  for (m = 0; m < 3*MAX_INTRAE_CONTRIBUTORS; m++) { 
    K_intracontrib.intraE_contributors_const[m] = (char) 5; 
  }

  kernelconstant_intra K_intra;
  for (m = 0; m < ATYPE_NUM; m++)	{ 
    K_intra.reqm_const[m] = (float) 0.1f; 
    K_intra.reqm_hbond_const[m] = (float) 0.1f;
    K_intra.atom1_types_reqm_const[m] = (unsigned int) 5;
    K_intra.atom2_types_reqm_const[m] = (unsigned int) 5;
  }
	for (m=0;m<MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES;m++)	{ 
    K_intra.VWpars_AC_const[m] = (float) 0.1f;
    K_intra.VWpars_BD_const[m] = (float) 0.1f;
  }
	for (m=0;m<MAX_NUM_OF_ATYPES;m++) {
    K_intra.dspars_S_const[m] = (float) 0.1f; 
    K_intra.dspars_V_const[m] = (float) 0.1f;
  }

	kernelconstant_rotlist K_rotlist;
  for (m = 0; m < MAX_NUM_OF_ROTATIONS; m++) { 
    K_rotlist.rotlist_const[m] = (int) 5; 
  }

	kernelconstant_conform K_conform;
  for (m = 0; m < MAX_NUM_OF_ATOMS; m++) { 
    K_conform.ref_coords_x_const[m]	= (float) 0.1f; 
    K_conform.ref_coords_y_const[m]	= (float) 0.1f; 
    K_conform.ref_coords_z_const[m] = (float) 0.1f;
  }
	for (m = 0; m < 3*MAX_NUM_OF_ROTBONDS; m++){ 
    K_conform.rotbonds_moving_vectors_const[m] = (float) 0.1f;
    K_conform.rotbonds_unit_vectors_const[m]  = (float) 0.1f;
  }
	for (m = 0; m < 4*MAX_NUM_OF_RUNS; m++) { 
    K_conform.ref_orientation_quats_const[m] = (float) 0.1f;
  }

  uint32_t NUM_MAX_ITERS = 2;

  #if defined (ENABLE_TRACE)
  ftrace_region_begin("OUTSIDE_LOOP_KERNEL1");
  #endif
  libkernel1( 
              23, // dockpars.num_of_atoms   
              6, // dockpars.num_of_atypes
              88, // dockpars.num_of_intraE_contributors
              81, // dockpars.gridsize_x
              81, // dockpars.gridsize_y
              81, // dockpars.gridsize_z
              6561, // g2
              531441, // g3
              0.37500, // dockpars.grid_spacing 
              105621835753504, // mem_dockpars_fgrids
              128, // dockpars.rotbondlist_length
              46.68814, // dockpars.coeff_elec 
              0.13220, // dockpars.coeff_desolv
              105621852759632, // mem_dockpars_conformations_current
              105621856599648, // mem_dockpars_energies_current
              105621860559696, // mem_dockpars_evals_of_new_entities
              150, // dockpars.pop_size
              0.01097, // dockpars.qasp 
              0.50000, // dockpars.smooth
              &K_interintra, //0x7ffe4b655370, // (char*)(&KerConst_interintra)  
              &K_intracontrib, //0x7ffe4b691e20, // (char*)(&KerConst_intracontrib)
              &K_intra, //0x7ffe4b655870, // (char*)(&KerConst_intra)   
              &K_rotlist, //0x7ffe4b65ade0, // (char*)(&KerConst_rotlist) 
              &K_conform, //0x7ffe4b656060, // (char*)(&KerConst_conform)a
              15000 // blocksPerGridForEachEntity
            );
  #if defined(ENABLE_TRACE)
  ftrace_region_end("OUTSIDE_LOOP_KERNEL1");
  #endif

  #if defined (ENABLE_TRACE)
  ftrace_region_begin("OUTSIDE_LOOP_KERNEL2");
  #endif
  libkernel2(
              150, // dockpars.pop_size 
              105621860559696, // mem_dockpars_evals_of_new_entities
              105621860619712, // mem_gpu_evals_of_runs
              100 // blocksPerGridForEachRun
            );
  #if defined(ENABLE_TRACE)
  ftrace_region_end("OUTSIDE_LOOP_KERNEL2");
  #endif

  #if defined (ENABLE_TRACE)
  ftrace_region_begin("KERNEL4");
  #endif
  for (uint32_t i=0; i<NUM_MAX_ITERS; i++) {
    libkernel4(
                23, // dockpars.num_of_atoms
                6, // dockpars.num_of_atypes
                88, // dockpars.num_of_intraE_contributors 
                81, // dockpars.gridsize_x
                81, // dockpars.gridsize_y
                81, // dockpars.gridsize_z
                6561, // g2
                531441, // g3 
                0.37500, // dockpars.grid_spacing
                105621835753504, // mem_dockpars_fgrids
                128, // dockpars.rotbondlist_length
                46.68814, // dockpars.coeff_elec
                0.13220, // dockpars.coeff_desolv
                105621852759632, // mem_dockpars_conformations_current
                105621856599648, // mem_dockpars_energies_current
                105621856659664, // mem_dockpars_conformations_next
                105621860499680, // mem_dockpars_energies_next
                105621860559696, // mem_dockpars_evals_of_new_entities
                105621860620128, // mem_dockpars_prng_states 
                150, // dockpars.pop_size
                8, // dockpars.num_of_genes
                0.60000, // dockpars.tournament_rate
                0.80000, // dockpars.crossover_rate
                0.02000, // dockpars.mutation_rate 
                16.00000, // dockpars.abs_max_dmov
                90.00000, // dockpars.abs_max_dang
                0.01097, // dockpars.qasp   
                0.50000, // dockpars.smooth 
                &K_interintra, //0x7ffe4b655370, // (char*)(&KerConst_interintra)  
                &K_intracontrib, //0x7ffe4b691e20, // (char*)(&KerConst_intracontrib)
                &K_intra, //0x7ffe4b655870, // (char*)(&KerConst_intra)   
                &K_rotlist, //0x7ffe4b65ade0, // (char*)(&KerConst_rotlist) 
                &K_conform, //0x7ffe4b656060, // (char*)(&KerConst_conform)a
                15000 // blocksPerGridForEachEntity
    );
  }
  #if defined(ENABLE_TRACE)
  ftrace_region_end("KERNEL4");
  #endif


  #if defined (ENABLE_TRACE)
  ftrace_region_begin("KERNEL3");
  #endif
  for (uint32_t i=0; i<NUM_MAX_ITERS; i++) {
    libkernel3(
                23, // dockpars.num_of_atoms
                6, // dockpars.num_of_atypes
                88, // dockpars.num_of_intraE_contributors
                81, // dockpars.gridsize_x
                81, // dockpars.gridsize_y
                81, // dockpars.gridsize_z
                6561, // g2
                531441, // g3
                0.37500, // dockpars.grid_spacing
                105621835753504, // mem_dockpars_fgrids
                128, // dockpars.rotbondlist_length
                46.68814, // dockpars.coeff_elec
                0.13220, // dockpars.coeff_desolv
                105621856659664, // mem_dockpars_conformations_next
                105621860499680, // mem_dockpars_energies_next
                105621860559696, // mem_dockpars_evals_of_new_entities
                105621860620128, // mem_dockpars_prng_states
                150, // dockpars.pop_size
                8, // dockpars.num_of_genes
                100.00000, // dockpars.lsearch_rate
                150, // dockpars.num_of_lsentities 
                0.01000, // dockpars.rho_lower_bound
                9.23760, // dockpars.base_dmov_mul_sqrt3
                129.90381, // dockpars.base_dang_mul_sqrt3
                4, // dockpars.cons_limit
                300, // dockpars.max_num_of_iters
                0.01097, // dockpars.qasp
                0.50000, // dockpars.smooth
                &K_interintra, //0x7ffe4b655370, // (char*)(&KerConst_interintra)  
                &K_intracontrib, //0x7ffe4b691e20, // (char*)(&KerConst_intracontrib)
                &K_intra, //0x7ffe4b655870, // (char*)(&KerConst_intra)   
                &K_rotlist, //0x7ffe4b65ade0, // (char*)(&KerConst_rotlist) 
                &K_conform, //0x7ffe4b656060, // (char*)(&KerConst_conform)a     
                15000 // blocksPerGridForEachLSEntity 
    );
  }
  #if defined(ENABLE_TRACE)
  ftrace_region_end("KERNEL3");
  #endif


  #if defined (ENABLE_TRACE)
  ftrace_region_begin("KERNEL2");
  #endif
  for (uint32_t i=0; i<NUM_MAX_ITERS; i++) {
    libkernel2(
                150, // dockpars.pop_size 
                105621860559696, // mem_dockpars_evals_of_new_entities
                105621860619712, // mem_gpu_evals_of_runs
                100 // blocksPerGridForEachRun
              );
  }
  #if defined(ENABLE_TRACE)
  ftrace_region_end("KERNEL2");
  #endif
*/
}
