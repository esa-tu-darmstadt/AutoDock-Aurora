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

typedef struct
{
	  float			    atom_charges_const	    	  [MAX_NUM_OF_ATOMS];
    int  			    atom_types_const  	     	  [MAX_NUM_OF_ATOMS];
	  int				    intraE_contributors_const   [3*MAX_INTRAE_CONTRIBUTORS];
    float 			  reqm_const 		     		      [ATYPE_NUM];
    float 			  reqm_hbond_const 	     	    [ATYPE_NUM];
    unsigned int 	atom1_types_reqm_const 	     [ATYPE_NUM];
    unsigned int	atom2_types_reqm_const 	     [ATYPE_NUM];
    float     		VWpars_AC_const              [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
    float     		VWpars_BD_const              [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
    float     		dspars_S_const               [MAX_NUM_OF_ATYPES];
    float     		dspars_V_const               [MAX_NUM_OF_ATYPES];
    int       		rotlist_const                [MAX_NUM_OF_ROTATIONS];
	  int       		subrotlist_1_const [MAX_NUM_OF_ROTATIONS];
	  int       		subrotlist_2_const [MAX_NUM_OF_ROTATIONS];
	  int       		subrotlist_3_const [MAX_NUM_OF_ROTATIONS];
	  int       		subrotlist_4_const [MAX_NUM_OF_ROTATIONS];
	  int       		subrotlist_5_const [MAX_NUM_OF_ROTATIONS];
	  int       		subrotlist_6_const [MAX_NUM_OF_ROTATIONS];
    int       		subrotlist_7_const [MAX_NUM_OF_ROTATIONS];
    int       		subrotlist_8_const [MAX_NUM_OF_ROTATIONS];
    int       		subrotlist_9_const [MAX_NUM_OF_ROTATIONS];
    float 			  ref_coords_x_const			 [MAX_NUM_OF_ATOMS];
	  float 			  ref_coords_y_const			 [MAX_NUM_OF_ATOMS];
	  float 			  ref_coords_z_const			 [MAX_NUM_OF_ATOMS];
    float			    rotbonds_moving_vectors_const[3*MAX_NUM_OF_ROTBONDS];
    float 			  rotbonds_unit_vectors_const  [3*MAX_NUM_OF_ROTBONDS];
	  float 			  ref_orientation_quats_const  [4*MAX_NUM_OF_RUNS];
} kernelconstant_static;

int main() {
  printf("This is main()\n");
  hello();


  // Using real VE VMA addresses for artifially invoking kernel2
  
  // Defining and initializing structs
  // as local variables -> will be stored in the stack!
  // This mimic what we did in our VEOffload design
  // as we passed some into the stack.
  // Initilization reference: host/src/calcenergy.cpp

  // init populations
  float cpu_init_populations[MAX_NUM_OF_RUNS][MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];

  for (int i = 0; i < MAX_NUM_OF_RUNS; i++) {
    for (int j = 0; j < MAX_POPSIZE; j++) {
      for (int k = 0; k < ACTUAL_GENOTYPE_LENGTH; k++) {
        cpu_init_populations[i][j][k] = 10.0f;
      }
    }
  }
  
  // prng
  unsigned int cpu_prng_seeds[MAX_NUM_OF_RUNS][MAX_POPSIZE];

  for (int i = 0; i < MAX_NUM_OF_RUNS; i++) {
    for (int j = 0; j < MAX_POPSIZE; j++) {
      cpu_prng_seeds[i][j] = 30;
    }
  }

  int m;
  kernelconstant_static KerConstStatic;

  for (m = 0; m < MAX_NUM_OF_ROTATIONS; m++) { 
    KerConstStatic.rotlist_const[m] = (int) 5;
    KerConstStatic.subrotlist_1_const[m] = (int) 5;
    KerConstStatic.subrotlist_2_const[m] = (int) 5;
    KerConstStatic.subrotlist_3_const[m] = (int) 5;
    KerConstStatic.subrotlist_4_const[m] = (int) 5;
    KerConstStatic.subrotlist_5_const[m] = (int) 5;
    KerConstStatic.subrotlist_6_const[m] = (int) 5;
    KerConstStatic.subrotlist_7_const[m] = (int) 5;
    KerConstStatic.subrotlist_8_const[m] = (int) 5;
    KerConstStatic.subrotlist_9_const[m] = (int) 5;
  }

  for (m = 0; m < MAX_NUM_OF_ATOMS; m++) { 
    KerConstStatic.ref_coords_x_const[m] = (float) 0.1f; 
    KerConstStatic.ref_coords_y_const[m] = (float) 0.1f; 
    KerConstStatic.ref_coords_z_const[m] = (float) 0.1f;
  }

	for (m = 0; m < 3*MAX_NUM_OF_ROTBONDS; m++){ 
    KerConstStatic.rotbonds_moving_vectors_const[m] = (float) 0.1f;
    KerConstStatic.rotbonds_unit_vectors_const[m]  = (float) 0.1f;
  }

	for (m = 0; m < 4*MAX_NUM_OF_RUNS; m++) { 
    KerConstStatic.ref_orientation_quats_const[m] = (float) 0.1f;
  }

  for (m = 0; m < MAX_NUM_OF_ATOMS; m++) { 
    KerConstStatic.atom_charges_const[m] = (float) 0.1f; 
    KerConstStatic.atom_types_const[m] = (char) 5;
  }

  for (m = 0; m < 3*MAX_INTRAE_CONTRIBUTORS; m++) { 
    KerConstStatic.intraE_contributors_const[m] = (char) 5; 
  }

  for (m = 0; m < ATYPE_NUM; m++)	{ 
    KerConstStatic.reqm_const[m] = (float) 0.1f; 
    KerConstStatic.reqm_hbond_const[m] = (float) 0.1f;
    KerConstStatic.atom1_types_reqm_const[m] = (unsigned int) 5;
    KerConstStatic.atom2_types_reqm_const[m] = (unsigned int) 5;
  }

	for (m=0;m<MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES;m++)	{ 
    KerConstStatic.VWpars_AC_const[m] = (float) 0.1f;
    KerConstStatic.VWpars_BD_const[m] = (float) 0.1f;
  }

	for (m=0;m<MAX_NUM_OF_ATYPES;m++) {
    KerConstStatic.dspars_S_const[m] = (float) 0.1f; 
    KerConstStatic.dspars_V_const[m] = (float) 0.1f;
  }

  #if defined (ENABLE_TRACE)
  ftrace_region_begin("OUTSIDE_LOOP_KERNEL_GA");
  #endif
  libkernel_ga(
              &cpu_init_populations[0][0][0],
              105621851785264, // mem_dockpars_conformations_current_Final
              105621851808080, // mem_dockpars_energies_current
              105621851808688, // mem_evals_performed
              105621851808720, // mem_gens_performed
              &cpu_prng_seeds[0][0],
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

              &KerConstStatic.rotlist_const[0],
              &KerConstStatic.subrotlist_1_const[0],
              &KerConstStatic.subrotlist_2_const[0],
              &KerConstStatic.subrotlist_3_const[0],
              &KerConstStatic.subrotlist_4_const[0],
              &KerConstStatic.subrotlist_5_const[0],
              &KerConstStatic.subrotlist_6_const[0],
              &KerConstStatic.subrotlist_7_const[0],
              &KerConstStatic.subrotlist_8_const[0],
              &KerConstStatic.subrotlist_9_const[0],
              &KerConstStatic.ref_coords_x_const[0],
              &KerConstStatic.ref_coords_y_const[0],
              &KerConstStatic.ref_coords_z_const[0],
              &KerConstStatic.rotbonds_moving_vectors_const[0],
              &KerConstStatic.rotbonds_unit_vectors_const[0],
              &KerConstStatic.ref_orientation_quats_const[0],
              31, // dockpars.rotbondlist_length
              23, // subrotlist_1_length
              8, // subrotlist_2_length
              0, // subrotlist_3_length
              0, // subrotlist_4_length
              0, // subrotlist_5_length
              0, // subrotlist_6_length
              0, // subrotlist_7_length
              0, // subrotlist_8_length
              0, // subrotlist_9_length

              &KerConstStatic.atom_charges_const[0],
              &KerConstStatic.atom_types_const[0],
              &KerConstStatic.intraE_contributors_const[0],
              &KerConstStatic.reqm_const[0],
              &KerConstStatic.reqm_hbond_const[0],
              &KerConstStatic.atom1_types_reqm_const[0],
              &KerConstStatic.atom2_types_reqm_const[0],
              &KerConstStatic.VWpars_AC_const[0],
              &KerConstStatic.VWpars_BD_const[0],
              &KerConstStatic.dspars_S_const[0],
              &KerConstStatic.dspars_V_const[0],
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
