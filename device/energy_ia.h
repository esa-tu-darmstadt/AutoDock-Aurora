#ifndef _INCLUDE_ENERGY_IA_H_
#define _INCLUDE_ENERGY_IA_H_

// --------------------------------------------------------------------------
// Calculates the intramolecular energy of a set of atomic ligand
// contributor-pairs.
// --------------------------------------------------------------------------
void energy_ia (
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
	 	float			DockConst_smooth,
		uint			DockConst_num_of_intraE_contributors,
		float			DockConst_grid_spacing,
		uint			DockConst_num_of_atypes,
		float			DockConst_coeff_elec,
		float			DockConst_qasp,
		float			DockConst_coeff_desolv,

        const	uint			DockConst_pop_size,
		float*			final_intraE,

		float	 	local_coords_x[][MAX_POPSIZE],
		float	 	local_coords_y[][MAX_POPSIZE],
		float 	 	local_coords_z[][MAX_POPSIZE]
                );
#endif
