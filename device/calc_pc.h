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

#ifndef _INCLUDE_CALC_PC_H_
#define _INCLUDE_CALC_PC_H_

// --------------------------------------------------------------------------
// Calculation of the pose of the ligand according to the input genotype.
// --------------------------------------------------------------------------
void calc_pc (
	const	int*	restrict	PC_rotlist,
	const	float*	restrict	PC_ref_coords_x,
	const	float*	restrict	PC_ref_coords_y,
	const	float*	restrict	PC_ref_coords_z,
	const	float*	restrict	PC_rotbonds_moving_vectors,
	const	float*	restrict	PC_rotbonds_unit_vectors,
	const	float*	restrict	PC_ref_orientation_quats,
	uint				DockConst_rotbondlist_length,
	uchar				DockConst_num_of_genes,
	uint				Host_RunId,

        const	uint			DockConst_pop_size,
		float		genotype[][MAX_POPSIZE],
		float	 	local_coords_x[][MAX_POPSIZE],
		float		local_coords_y[][MAX_POPSIZE],
		float		local_coords_z[][MAX_POPSIZE]
              );
#endif
