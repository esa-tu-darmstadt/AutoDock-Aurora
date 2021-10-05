/*

AutoDock-Aurora, a vectorized implementation of AutoDock 4.2 for the NEC SX-Aurora TSUBASA
Copyright (C) 2021 TU Darmstadt, Embedded Systems and Applications Group, Germany. All rights reserved.
For some of the code, Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

AutoDock is a Trade Mark of the Scripps Research Institute.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef _INCLUDE_LGA_H_
#define _INCLUDE_LGA_H_

// --------------------------------------------------------------------------
// Lamarckian Genetic-Algorithm (GA): GA + LS (Local Search)
// --------------------------------------------------------------------------
void lga (
	const 	float*		PopulationCurrentInitial,
			float*		PopulationCurrentFinal,
			float*		EnergyCurrent,
			uint*		Evals_performed,
			uint*		Gens_performed,
			uint		DockConst_pop_size,
			uint		DockConst_num_of_energy_evals,
			uint		DockConst_num_of_generations,
			float		DockConst_tournament_rate,
			float		DockConst_mutation_rate,
			float		DockConst_abs_max_dmov,
			float		DockConst_abs_max_dang,
			float		Host_two_absmaxdmov,
			float		Host_two_absmaxdang,
			float		DockConst_crossover_rate,
			uchar		DockConst_num_of_genes,
	// PC
	const	int* 		PC_rotlist,
	const	float*		PC_ref_coords_x, // TODO: merge them into a single one?
	const	float*		PC_ref_coords_y,
	const	float*		PC_ref_coords_z,
	const	float*		PC_rotbonds_moving_vectors,
	const	float*		PC_rotbonds_unit_vectors,
	const	float*		PC_ref_orientation_quats,
			uint		DockConst_rotbondlist_length,
	// IA
	const 	float*		IA_IE_atom_charges,
	const	int*		IA_IE_atom_types,
	const	int*		IA_intraE_contributors,
	const	float*		IA_reqm,
	const	float*		IA_reqm_hbond,
	const	uint*		IA_atom1_types_reqm,
	const	uint*		IA_atom2_types_reqm,
	const	float*		IA_VWpars_AC,
	const	float*		IA_VWpars_BD,
	const	float*		IA_dspars_S,
	const	float*		IA_dspars_V,
			float		DockConst_smooth,
			uint		DockConst_num_of_intraE_contributors,
			float		DockConst_grid_spacing,
			uint		DockConst_num_of_atypes,
			float		DockConst_coeff_elec,
			float		DockConst_qasp,
			float		DockConst_coeff_desolv,
	// IE
	const	float*		Fgrids,
			uchar		DockConst_xsz,
			uchar		DockConst_ysz,
			uchar		DockConst_zsz,
			uchar		DockConst_num_of_atoms,
			float		DockConst_gridsize_x_minus1,
			float		DockConst_gridsize_y_minus1,
			float		DockConst_gridsize_z_minus1,
			uint		Host_mul_tmp2,
			uint		Host_mul_tmp3,
	// LS
			uchar		lsmet,
	// LS-SW
			ushort		DockConst_max_num_of_iters,
			float		DockConst_rho_lower_bound,
			float		DockConst_base_dmov_mul_sqrt3,
			float		DockConst_base_dang_mul_sqrt3,
			uchar		DockConst_cons_limit,
	// LS-AD
	const 	int*		GRAD_rotbonds,
	const 	int*		GRAD_rotbonds_atoms,
	const 	int*		GRAD_num_rotating_atoms_per_rotbond,
	const 	float*		GRAD_angle,
	const 	float*		GRAD_dependence_on_theta,
	const 	float*		GRAD_dependence_on_rotangle,
	// Values changing every LGA run
			uint		Host_RunId,
			uint 	    Host_Offset_Pop,
			uint	    Host_Offset_Ene
    );

#endif
