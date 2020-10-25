#ifndef CALCENERGY_H_
#define CALCENERGY_H_

#include "ext_headers.h"
#include "miscellaneous.h"
#include "processligand.h"
#include "getparameters.h"

// Members of this struct are passed to the FPGA OpenCL kernels as inputs.
// Struct members are parameters related to the ligand, the grid
// and the genetic algorithm, or they are pointers of FPGA memory areas
// used for storing different data such as the current
// and the next population genotypes and energies, the grids,
// the evaluation counters and the random number generator states.
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
	float  			grid_spacing;
	unsigned int    rotbondlist_length;
	float  			coeff_elec;
	float  			coeff_desolv;
	unsigned int 	num_of_energy_evals;
	unsigned int 	num_of_generations;
	unsigned int    pop_size;
	unsigned char   num_of_genes;
	float  			tournament_rate;
	float  			crossover_rate;
	float  			mutation_rate;
	float  			abs_max_dmov;
	float  			abs_max_dang;
	float  			lsearch_rate;
	float 			smooth;
	unsigned int 	num_of_lsentities;
	float  			rho_lower_bound;
	float  			base_dmov_mul_sqrt3;
	float  			base_dang_mul_sqrt3;
	unsigned int 	cons_limit;
	unsigned int 	max_num_of_iters;
	float  			qasp;
} Dockparameters;

typedef struct
{
	float			atom_charges_const	    	 [MAX_NUM_OF_ATOMS];
    char  			atom_types_const  	     	 [MAX_NUM_OF_ATOMS];
	char			intraE_contributors_const    [3*MAX_INTRAE_CONTRIBUTORS];
    float 			reqm_const 		     		 [ATYPE_NUM];
    float 			reqm_hbond_const 	     	 [ATYPE_NUM];
    unsigned int 	atom1_types_reqm_const 	     [ATYPE_NUM];
    unsigned int	atom2_types_reqm_const 	     [ATYPE_NUM];
    float     		VWpars_AC_const              [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
    float     		VWpars_BD_const              [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
    float     		dspars_S_const               [MAX_NUM_OF_ATYPES];
    float     		dspars_V_const               [MAX_NUM_OF_ATYPES];
    int       		rotlist_const                [MAX_NUM_OF_ROTATIONS];
    float 			ref_coords_x_const			 [MAX_NUM_OF_ATOMS];
	float 			ref_coords_y_const			 [MAX_NUM_OF_ATOMS];
	float 			ref_coords_z_const			 [MAX_NUM_OF_ATOMS];
    float			rotbonds_moving_vectors_const[3*MAX_NUM_OF_ROTBONDS;
    float 			rotbonds_unit_vectors_const  [3*MAX_NUM_OF_ROTBONDS];
	float 			ref_orientation_quats_const  [MAX_NUM_OF_RUNS][4];
} kernelconstant_static;

int prepare_conststatic_fields_for_fpga(Liganddata* 	       myligand_reference,
				 	Dockpars*   	       mypars,
				 	float*      	       cpu_ref_ori_angles,
				 	kernelconstant_static* KerConstStatic);

void make_reqrot_ordering(char number_of_req_rotations[MAX_NUM_OF_ATOMS],
			  char atom_id_of_numrots[MAX_NUM_OF_ATOMS],
		          int  num_of_atoms);

int gen_rotlist(Liganddata* myligand,
		int         rotlist[MAX_NUM_OF_ROTATIONS]);

#endif /* CALCENERGY_H_ */






