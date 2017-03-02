/*
 * (C) 2013. Evopro Innovation Kft.
 *
 * getparameters.h
 *
 *  Created on: 2008.10.22.
 *      Author: pechan.imre
 */

#ifndef GETPARAMETERS_H_
#define GETPARAMETERS_H_

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "defines.h"
#include "processligand.h"
#include "processgrid.h"
#include "miscallenous.h"

/*
typedef struct
{
	double AD4_coeff_vdW;
	double AD4_coeff_hb;
	double scaled_AD4_coeff_elec;
	double AD4_coeff_desolv;
	double AD4_coeff_tors;
} AD4_free_energy_coeffs;
*/

typedef struct { 	
	float AD4_coeff_vdW; 	
	float AD4_coeff_hb; 	
	float scaled_AD4_coeff_elec; 	
	float AD4_coeff_desolv; 	
	float AD4_coeff_tors; 
} AD4_free_energy_coeffs;

/*
typedef struct
//Struct which contains the docking parameters (partly parameters for fpga)
{
	unsigned long num_of_genes;
	unsigned long num_of_atoms;
	unsigned long grid_size_xyz [3];
	unsigned long sq_spacing;
	unsigned long num_of_energy_evals;
	unsigned long num_of_generations;
	unsigned long seed;
	unsigned long dmov_mask;
	unsigned long dang_mask;
	unsigned long mutation_rate;
	unsigned long crossover_rate;
	double lsearch_rate;
	unsigned long num_of_ls;
	unsigned long tournament_rate;
	unsigned long rho_lower_bound;
	unsigned long base_dmov_mul_sqrt3;
	unsigned long base_dang_mul_sqrt3;
	unsigned long cons_limit;
	unsigned long max_num_of_iters;
	unsigned long pop_size;
	unsigned long submodules_enable;
	unsigned long size_of_cube;
	unsigned long size_of_plane;
	unsigned long size_of_line;
	unsigned long eldesgrid_id;
	unsigned long adr0;
	unsigned long adr1;
	unsigned long adr2;
	unsigned long adr3;
	unsigned long adr4;
	unsigned long adr5;
	char seed_gen_or_loadfile;
	char initpop_gen_or_loadfile;
	unsigned char gen_pdbs;
	char fldfile [128];
	char ligandfile [128];
	double ref_ori_angles [3];
	int num_of_runs;
	char reflig_en_reqired;
	unsigned long evals_performed;
	unsigned long generations_used;
	char unbound_model;
	AD4_free_energy_coeffs coeffs;
	char handle_symmetry;
	char gen_finalpop;
	char gen_best;
	char resname [128];
	double qasp;
	double rmsd_tolerance;
} Dockpars;
*/

typedef struct 
//Struct which contains the docking parameters (partly parameters for fpga) 
{ 	unsigned int num_of_genes; 	
	unsigned int num_of_atoms; 	
	unsigned int grid_size_xyz [3]; 	
	unsigned int sq_spacing; 	
	unsigned int num_of_energy_evals; 	
	unsigned int num_of_generations; 	
	unsigned int seed; 	
	unsigned int dmov_mask; 	
	unsigned int dang_mask; 	
	unsigned int mutation_rate; 	
	unsigned int crossover_rate; 	
	float lsearch_rate; 	
	unsigned int num_of_ls; 	
	unsigned int tournament_rate; 	
	unsigned int rho_lower_bound; 	
	unsigned int base_dmov_mul_sqrt3; 	
	unsigned int base_dang_mul_sqrt3; 	
	unsigned int cons_limit; 	
	unsigned int max_num_of_iters; 	
	unsigned int pop_size; 	
	unsigned int submodules_enable; 	
	unsigned int size_of_cube; 	
	unsigned int size_of_plane; 	
	unsigned int size_of_line; 	
	unsigned int eldesgrid_id; 	
	unsigned int adr0; 	
	unsigned int adr1; 	
	unsigned int adr2; 	
	unsigned int adr3; 	
	unsigned int adr4; 	
	unsigned int adr5; 	
	char seed_gen_or_loadfile; 	
	char initpop_gen_or_loadfile; 	
	unsigned char gen_pdbs; 	
	char fldfile [128]; 	
	char ligandfile [128]; 	
	float ref_ori_angles [3]; 	
	int num_of_runs; 	
	char reflig_en_reqired; 	
	unsigned int evals_performed; 	
	unsigned int generations_used; 	
	char unbound_model; 	
	AD4_free_energy_coeffs coeffs; 	
	char handle_symmetry; 	
	char gen_finalpop; 	
	char gen_best; 	
	char resname [128]; 	
	float qasp; 	
	float rmsd_tolerance; 
} Dockpars;

int get_filenames_and_ADcoeffs(const int*, char**, Dockpars*);

void get_recligpars(const Liganddata*, const Gridinfo*, Dockpars*);

/*
void get_commandpars(const int*, char**, double*, Dockpars*);
*/
void get_commandpars(const int*, 		    
		     char**, 		     
		     float*, 		     
		     Dockpars*);



/*
void get_seeds_and_initpop(Dockpars*, double [][40], const Liganddata*, const Gridinfo*);
*/
void get_seeds_and_initpop(Dockpars*, 			   
			   float [][40], 			   
			   const Liganddata*, 			   
			   const Gridinfo*);


void get_ref_orientation(Liganddata*, const Dockpars*);

#endif /* GETPARAMETERS_H_ */


