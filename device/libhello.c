#include <stdio.h>
#include <stdint.h>
#include <ftrace.h>

#include "defines.h"
#include "libkernel_ga.c"


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

*/


  #if defined (ENABLE_TRACE)
  ftrace_region_begin("OUTSIDE_LOOP_KERNEL_GA");
  #endif
  libkernel_ga(
              105621851762448, // mem_dockpars_conformations_current_Initial
              105621851785264, // mem_dockpars_conformations_current_Final
              105621851808080, // mem_dockpars_energies_current
              105621851808688, // mem_evals_performed
              105621851808720, // mem_gens_performed
              105621851808752, // mem_prng_states
              150, // dockpars.pop_size
              2048000, // dockpars.num_of_energy_evals
              99999, // dockpars.num_of_generations
              0.6, // dockpars.tournament_rate
              0.02, // dockpars.mutation_rate
              16, // dockpars.abs_max_dmov
              90, // dockpars.abs_max_dang
              32, // two_absmaxdmov
              180, // two_absmaxdang
              0.8, // dockpars.crossover_rate
              150, // dockpars.num_of_lsentities
              8, // dockpars.num_of_genes

              105621851809360, // mem_pc_rotlist_const
              105621851825760, // mem_pc_ref_coords_x_const
              105621851826800, // mem_pc_ref_coords_y_const
              105621851827840, // mem_pc_ref_coords_z_const
              105621851828880, // mem_pc_rotbonds_moving_vectors_const
              105621851829280, // mem_pc_rotbonds_unit_vectors_const
              105621851829680, // mem_pc_ref_orientation_quats_const
              31, // dockpars.rotbondlist_length

              105621868837424, // mem_ia_ie_atom_charges_const
              105621868838464, // mem_ia_ie_atom_types_const
              105621868838736, // mem_ia_contributors_const
              105621868936288, // mem_ia_reqm_const
              105621868936384, // mem_ia_reqm_hbond_const
              105621868936480, // mem_ia_atom1_types_reqm_const
              105621868936576, // mem_ia_atom2_types_reqm_const
              105621868936672, // mem_ia_VWpars_AC_const
              105621868937472, // mem_ia_VWpars_BD_const
              105621868938272, // mem_ia_dspars_S_const
              105621868938336, // mem_ia_dspars_V_const
              0.5, // dockpars.smooth
              88, // dockpars.num_of_intraE_contributors
              0.375, // dockpars.grid_spacing
              6, // dockpars.num_of_atypes
              46.6881, // dockpars.coeff_elec
              0.01097, // dockpars.qasp
              0.1322, // dockpars.coeff_desolv

              105621851831296, // mem_dockpars_fgrids
              81, // dockpars.g1
              6561, // dockpars.g2
              531441, // dockpars.g3
              23, // dockpars.num_of_atoms
              80, // fgridsizex_minus1
              80, // fgridsizey_minus1
              80, // fgridsizez_minus1
              3188646, // mul_tmp2
              3720087, // mul_tmp3

              300, // Host_max_num_of_iters
              0.01, // dockpars.rho_lower_bound
              9.2376, // dockpars.base_dmov_mul_sqrt3
              129.904, // dockpars.base_dang_mul_sqrt3
              4, // Host_cons_limit
              1 // mypars->num_of_runs
            );
  #if defined(ENABLE_TRACE)
  ftrace_region_end("OUTSIDE_LOOP_KERNEL_GA");
  #endif

}
