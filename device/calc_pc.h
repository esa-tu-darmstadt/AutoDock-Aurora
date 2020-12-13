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
