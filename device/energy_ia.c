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

#include "auxiliary.h"

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
	 	float				DockConst_smooth,
		uint				DockConst_num_of_intraE_contributors,
		float				DockConst_grid_spacing,
		uint				DockConst_num_of_atypes,
		float				DockConst_coeff_elec,
		float				DockConst_qasp,
		float				DockConst_coeff_desolv,

        const	uint				DockConst_pop_size,
		float*				final_intraE,

		float	 	local_coords_x[][MAX_POPSIZE],
		float	 	local_coords_y[][MAX_POPSIZE],
		float 	 	local_coords_z[][MAX_POPSIZE]
)
{
#if defined (PRINT_ALL_IA) 
	printf("\n");
	printf("Starting <intra-molecular calculation> ... \n");
	printf("\n");
	printf("\t%-40s %f\n", "DockConst_smooth: ",	DockConst_smooth);
	printf("\t%-40s %u\n", "DockConst_num_of_intraE_contributors: ",	DockConst_num_of_intraE_contributors);
	printf("\t%-40s %f\n", "DockConst_grid_spacing: ",	DockConst_grid_spacing);
	printf("\t%-40s %u\n", "DockConst_num_of_atypes: ",	DockConst_num_of_atypes);
	printf("\t%-40s %f\n", "DockConst_coeff_elec: ",	DockConst_coeff_elec);
	printf("\t%-40s %f\n", "DockConst_qasp: ",	DockConst_qasp);
	printf("\t%-40s %f\n", "DockConst_coeff_desolv: ",	DockConst_coeff_desolv);
#endif

#if defined (ENABLE_TRACE)
	ftrace_region_begin("IA_MAIN_LOOP");
#endif
	for (int j = 0; j < DockConst_pop_size; j++)
	{
		final_intraE[j] = 0.0f;      
	}

	// For each intramolecular atom contributor pair
	for (uint contributor_counter = 0; contributor_counter < DockConst_num_of_intraE_contributors; contributor_counter++)
	{
		int atom1_id = IA_intraE_contributors[3*contributor_counter];
		int atom2_id = IA_intraE_contributors[3*contributor_counter + 1];

		// Getting types ids
		int atom1_typeid = IA_IE_atom_types[atom1_id];
		int atom2_typeid = IA_IE_atom_types[atom2_id];

		// Getting optimum pair distance (opt_distance) from reqm and reqm_hbond
		// reqm: equilibrium internuclear separation 
		//       (sum of the vdW radii of two like atoms (A)) in the case of vdW
		// reqm_hbond: equilibrium internuclear separation
		// 	 (sum of the vdW radii of two like atoms (A)) in the case of hbond 
		float opt_distance;
		float delta_distance = 0.5f*DockConst_smooth;

		int atom1_type_vdw_hb = IA_atom1_types_reqm[atom1_typeid];
		int atom2_type_vdw_hb = IA_atom2_types_reqm[atom2_typeid];
			
		if (IA_intraE_contributors[3*contributor_counter + 2] == 1)	// H-bond
		{
			opt_distance = IA_reqm_hbond[atom1_type_vdw_hb] + IA_reqm_hbond[atom2_type_vdw_hb];
		}
		else	// Van der Waals
		{
			opt_distance = 0.5f*(IA_reqm[atom1_type_vdw_hb] + IA_reqm[atom2_type_vdw_hb]);
		}

		int is_H_bond = IA_intraE_contributors[3*contributor_counter + 2];
		float vdW_const1 = IA_VWpars_AC[atom1_typeid * DockConst_num_of_atypes + atom2_typeid];
		float vdW_const2 = IA_VWpars_BD[atom1_typeid * DockConst_num_of_atypes + atom2_typeid];

		float desolv_const = ((IA_dspars_S[atom1_typeid] + DockConst_qasp * esa_fabs(IA_IE_atom_charges[atom1_id])) * IA_dspars_V[atom2_typeid] +
				      (IA_dspars_S[atom2_typeid] + DockConst_qasp * esa_fabs(IA_IE_atom_charges[atom2_id])) * IA_dspars_V[atom1_typeid]) *
			DockConst_coeff_desolv;
		float elec_const = DockConst_coeff_elec * IA_IE_atom_charges[atom1_id] * IA_IE_atom_charges[atom2_id];

#pragma _NEC packed_vector
#pragma _NEC vovertake
                //#pragma _NEC advance_gather   # this directive is dangerous here!
#pragma _NEC gather_reorder
#pragma omp simd simdlen(512)
		for (int j = 0; j < DockConst_pop_size; j++)
		{
			float subx = local_coords_x[atom1_id][j] - local_coords_x[atom2_id][j];
			float suby = local_coords_y[atom1_id][j] - local_coords_y[atom2_id][j];
			float subz = local_coords_z[atom1_id][j] - local_coords_z[atom2_id][j];

			// Using built-int sqrt() is some seconds slower
			float atomic_distance = esa_sqrt(subx*subx + suby*suby + subz*subz)*DockConst_grid_spacing;

#if defined (PRINT_ALL)
			printf("\nContrib %u: pop %3d, atoms %u and %u, distance: %f\n", contributor_counter, j, (atom1_id+1), (atom2_id+1), atomic_distance);
#endif
		
			// But only if the distance is less than distance cutoff value and 20.48A (because of the tables)
			//if ((distance_leo < 8.0f) && (distance_leo < 20.48f))
			float partialE1 = 0.0f;
			float partialE2 = 0.0f;
			float partialE3 = 0.0f;
			float partialE4 = 0.0f;

			// Getting smoothed distance
			// smoothed_distance = function(atomic_distance, opt_distance)
			float smoothed_distance;

			if (atomic_distance <= (opt_distance - delta_distance)) {
				smoothed_distance = atomic_distance + delta_distance;
			}
			else if (atomic_distance < (opt_distance + delta_distance)) {
				smoothed_distance = opt_distance;
			}
			else { // else if (atomic_distance >= (opt_distance + delta_distance))
				smoothed_distance = atomic_distance - delta_distance;
			}

			float distance_pow_2  = atomic_distance*atomic_distance; 

			float inverse_smoothed_distance_pow_2  = 1.0f / (smoothed_distance*smoothed_distance);
			float inverse_smoothed_distance_pow_4  = inverse_smoothed_distance_pow_2 * inverse_smoothed_distance_pow_2;
			float inverse_smoothed_distance_pow_6  = inverse_smoothed_distance_pow_4 * inverse_smoothed_distance_pow_2;

			// Calculating energy contributions
			// Cuttoff1: internuclear-distance at 8A only for vdw and hbond.
			if (atomic_distance < 8.0f) 
			{
				// Calculating van der Waals / hydrogen bond term

				if (is_H_bond == 1) {	// H-bond
					float inverse_smoothed_distance_pow_10 = inverse_smoothed_distance_pow_6 * inverse_smoothed_distance_pow_4;
					float inverse_smoothed_distance_pow_12 = inverse_smoothed_distance_pow_6 * inverse_smoothed_distance_pow_6;
					partialE1 = vdW_const1 * inverse_smoothed_distance_pow_12;
					partialE2 = vdW_const2 * inverse_smoothed_distance_pow_10;
				}
				else {	// Van der Waals
					float inverse_smoothed_distance_pow_12 = inverse_smoothed_distance_pow_6 * inverse_smoothed_distance_pow_6;
					partialE1 = vdW_const1 * inverse_smoothed_distance_pow_12;
					partialE2 = vdW_const2 * inverse_smoothed_distance_pow_6;
				}
			} // if cuttoff1 - internuclear-distance at 8A

			// Calculating energy contributions
			// Cuttoff2: internuclear-distance at 20.48A only for el and sol.
			if (atomic_distance < 20.48f)
			{
				float term_partialE3 = atomic_distance *
					(DIEL_A + (DIEL_B / (1.0f + DIEL_K * esa_expf0(DIEL_B_TIMES_H * atomic_distance))));
				float term_inv_partialE3 = (1.0f / term_partialE3);

				// Calculating electrostatic term
				partialE3 = elec_const * term_inv_partialE3;
              
				// Calculating desolvation term
				partialE4 = desolv_const * esa_expf0(-0.0386f * distance_pow_2);
			} // if cuttoff2 - internuclear-distance at 20.48A
	
			final_intraE[j] += partialE1 - partialE2 + partialE3 + partialE4;
	
		} // End j Loop (over individuals)

	} // End of contributor_counter for-loop

#if defined (ENABLE_TRACE)
	ftrace_region_end("IA_MAIN_LOOP");
#endif

	// --------------------------------------------------------------

#if defined (PRINT_ALL_IA) 
	printf("\n");
	printf("Finishing <intra-molecular calculation>\n");
	printf("\n");
#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
