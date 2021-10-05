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

#ifndef _INCLUDE_PERFORM_LS_H_
#define _INCLUDE_PERFORM_LS_H_

#include "rand_gen.h"

void perform_ls(
			ushort				DockConst_max_num_of_iters,
			float				DockConst_rho_lower_bound,
			float				DockConst_base_dmov_mul_sqrt3,
			uchar				DockConst_num_of_genes,
			float				DockConst_base_dang_mul_sqrt3,
			uchar				DockConst_cons_limit,
			uint				pop_size,
			float				in_out_genotype[][MAX_POPSIZE],
			float*	restrict	in_out_energy,
			uint*	restrict	out_eval,
	// PC
	const	int*	restrict	PC_rotlist,
	const	float*	restrict	PC_ref_coords_x,
	const	float*	restrict	PC_ref_coords_y,
	const	float*	restrict	PC_ref_coords_z,
	const	float*	restrict	PC_rotbonds_moving_vectors,
	const	float*	restrict	PC_rotbonds_unit_vectors,
	const	float*	restrict	PC_ref_orientation_quats,
			uint				DockConst_rotbondlist_length,
			uint				Host_RunId,
	// IA
	const	float*	restrict	IA_IE_atom_charges,
	const	int*	restrict	IA_IE_atom_types,
	const	int*	restrict	IA_intraE_contributors,
	const	float*	restrict	IA_reqm,
	const	float*	restrict	IA_reqm_hbond,
	const	uint*	restrict	IA_atom1_types_reqm,
	const	uint*	restrict	IA_atom2_types_reqm,
	const	float*	restrict	IA_VWpars_AC,
	const	float*	restrict	IA_VWpars_BD,
	const	float*	restrict	IA_dspars_S,
	const	float*	restrict	IA_dspars_V,
			float				DockConst_smooth,
			uint				DockConst_num_of_intraE_contributors,
			float				DockConst_grid_spacing,
			uint				DockConst_num_of_atypes,
			float				DockConst_coeff_elec,
			float				DockConst_qasp,
			float				DockConst_coeff_desolv,
	// IE
	const	float*	restrict	IE_Fgrids,
			uchar				DockConst_xsz,
			uchar				DockConst_ysz,
			uchar				DockConst_zsz,
			uchar				DockConst_num_of_atoms,
			float				DockConst_gridsize_x_minus1,
			float				DockConst_gridsize_y_minus1,
			float				DockConst_gridsize_z_minus1,
			uint				Host_mul_tmp2,
			uint				Host_mul_tmp3
                );
#endif
