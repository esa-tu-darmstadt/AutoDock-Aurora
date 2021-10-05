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

#ifndef _INCLUDE_ENERGY_IE_H_
#define _INCLUDE_ENERGY_IE_H_

// --------------------------------------------------------------------------
// Calculates the intermolecular energy of a ligand given by
// ligand xyz-positions, and a receptor represented as a grid. 
// The grid point values must be stored at the location which starts at GlobFgrids. 
// If an atom is remains outside the grid, 
// a very high value will be added to the current energy as a penalty. 
// --------------------------------------------------------------------------
void energy_ie (
	const	float*	restrict	IE_Fgrids,
	const	float*	restrict	IA_IE_atom_charges,
	const	int*	restrict	IA_IE_atom_types,
			uchar				DockConst_xsz,
			uchar				DockConst_ysz,
			uchar				DockConst_zsz,
			uchar				DockConst_num_of_atoms,
			float				DockConst_gridsize_x_minus1,
			float				DockConst_gridsize_y_minus1,
			float				DockConst_gridsize_z_minus1,
			uint				Host_mul_tmp2,
			uint				Host_mul_tmp3,

	        const	uint				DockConst_pop_size,
			float*				final_interE,

			float		local_coords_x[][MAX_POPSIZE],
			float		local_coords_y[][MAX_POPSIZE],
			float		local_coords_z[][MAX_POPSIZE]
                );
#endif
