#ifndef DEVICE_ARGS_H_
#define DEVICE_ARGS_H

#ifndef SHORT_TYPE_NAMES_
#define SHORT_TYPE_NAMES_
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int  uint;
#endif

struct device_args {
	int  ve_proc_id;
	int  ve_num_procs;
	const float *PopulationCurrentInitial;
	float *PopulationCurrentFinal;
	float *EnergyCurrent;
	uint  *Evals_performed;
	uint  *Gens_performed;
	uint  *dockpars_prng_states;
	uint  DockConst_pop_size;
	uint  DockConst_num_of_energy_evals;
	uint  DockConst_num_of_generations;
	float DockConst_tournament_rate;
	float DockConst_mutation_rate;
	float DockConst_abs_max_dmov;
	float DockConst_abs_max_dang;
	float Host_two_absmaxdmov;
	float Host_two_absmaxdang;
	float DockConst_crossover_rate;
	uint  DockConst_num_of_lsentities;
	uchar DockConst_num_of_genes;
	// PC
	const int   *PC_rotlist;
	const float *PC_ref_coords_x;// TODO: merge them into a single one?
	const float *PC_ref_coords_y;
	const float *PC_ref_coords_z;
	const float *PC_rotbonds_moving_vectors;
	const float *PC_rotbonds_unit_vectors;
	const float *PC_ref_orientation_quats;
          uint  DockConst_rotbondlist_length;
	// IA
	const float *IA_IE_atom_charges;
	const int   *IA_IE_atom_types;
	const int   *IA_intraE_contributors;
	const float *IA_reqm;
	const float *IA_reqm_hbond;
	const uint  *IA_atom1_types_reqm;
	const uint  *IA_atom2_types_reqm;
	const float *IA_VWpars_AC;
	const float *IA_VWpars_BD;
	const float *IA_dspars_S;
	const float *IA_dspars_V;
	float DockConst_smooth;
	uint  DockConst_num_of_intraE_contributors;
	float DockConst_grid_spacing;
	uint  DockConst_num_of_atypes;
	float DockConst_coeff_elec;
	float DockConst_qasp;
	float DockConst_coeff_desolv;
	// IE
	const float *Fgrids;
	uchar DockConst_xsz;
	uchar DockConst_ysz;
	uchar DockConst_zsz;
	uchar DockConst_num_of_atoms;
	float DockConst_gridsize_x_minus1;
	float DockConst_gridsize_y_minus1;
	float DockConst_gridsize_z_minus1;
	uint  Host_mul_tmp2;
	uint  Host_mul_tmp3;
	// Solis-Wets (LS-SW)
	ushort DockConst_max_num_of_iters;
	float  DockConst_rho_lower_bound;
	float  DockConst_base_dmov_mul_sqrt3;
	float  DockConst_base_dang_mul_sqrt3;
	uchar  DockConst_cons_limit;
	// ADADELTA (LS-AD)

	// Values changing every LGA run
	uint   Host_num_of_runs;
};

#endif // DEVICE_ARGS_H_
