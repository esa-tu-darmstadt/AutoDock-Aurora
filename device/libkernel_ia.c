#include "defines.h"
#include "math.h"
#include "auxiliary.c"

// --------------------------------------------------------------------------
// IntraE calculates the intramolecular energy of a set of atomic ligand 
// contributor-pairs.
// Originally from: processligand.c
// --------------------------------------------------------------------------
void libkernel_ia (
	const   float* restrict KerConstStatic_atom_charges_const,
 	const	char*  restrict KerConstStatic_atom_types_const,

	const 	char*  restrict KerConstStatic_intraE_contributors_const,

			float 				 DockConst_smooth,
	const   float* restrict KerConstStatic_reqm,
	const   float* restrict KerConstStatic_reqm_hbond,    
	const   uint*  restrict KerConstStatic_atom1_types_reqm,
	const   uint*  restrict KerConstStatic_atom2_types_reqm,  

	const   float* restrict KerConstStatic_VWpars_AC_const,
	const   float* restrict KerConstStatic_VWpars_BD_const,
	const   float* restrict KerConstStatic_dspars_S_const,
 	const   float* restrict KerConstStatic_dspars_V_const,

		   	uint			DockConst_num_of_intraE_contributors,
		  	float           DockConst_grid_spacing,
			uchar           DockConst_num_of_atypes,
			float           DockConst_coeff_elec,
			float           DockConst_qasp,
			float           DockConst_coeff_desolv,
			float* 			final_intraE,

			float* 			restrict local_coords_x,
			float* 			restrict local_coords_y,
			float* 			restrict local_coords_z
)
{
	#if defined (DEBUG_ACTIVE_KERNEL)
	if (active == 0) {printf("	%-20s: %s\n", "Krnl_IntraE", "must be disabled");}
	#endif

	float intraE = 0.0f;

	// For each intramolecular atom contributor pair

	// LOOP_FOR_INTRAE_MAIN
	for (ushort contributor_counter=0; contributor_counter<DockConst_num_of_intraE_contributors; contributor_counter++) {

		char atom1_id = KerConstStatic_intraE_contributors_const[3*contributor_counter];
		char atom2_id = KerConstStatic_intraE_contributors_const[3*contributor_counter + 1];

		float subx = local_coords_x[atom1_id] - local_coords_x[atom2_id];
		float suby = local_coords_y[atom1_id] - local_coords_y[atom2_id];
		float subz = local_coords_z[atom1_id] - local_coords_z[atom2_id];

		//atomic_distance = sqrt(subx*subx + suby*suby + subz*subz)*DockConst_grid_spacing;
		float atomic_distance = sqrt_custom(subx*subx + suby*suby + subz*subz)*DockConst_grid_spacing;

/*
		if (atomic_distance < 1.0f) {
			#if defined (DEBUG_KRNL_INTRAE)
			printf("\n\nToo low distance (%f) between atoms %u and %u\n", atomic_distance, atom1_id, atom2_id);
			#endif
			//return HIGHEST_ENERGY;	// Returning maximal value
			atomic_distance = 1.0f;
		}
*/

		#if defined (DEBUG_KRNL_INTRAE)
		printf("\n\nCalculating energy contribution of atoms %u and %u\n", atom1_id+1, atom2_id+1);
		printf("Distance: %f\n", atomic_distance);
		#endif

/*
		float partialE1;
		float partialE2;
		float partialE3;
		float partialE4;
*/
		// But only if the distance is less than distance cutoff value and 20.48A (because of the tables)
		//if ((distance_leo < 8.0f) && (distance_leo < 20.48f))
		float partialE1 = 0.0f;
		float partialE2 = 0.0f;
		float partialE3 = 0.0f;
		float partialE4 = 0.0f;

 		// Getting types ids
		char atom1_typeid = KerConstStatic_atom_types_const [atom1_id];
		char atom2_typeid = KerConstStatic_atom_types_const [atom2_id];

		// Getting optimum pair distance (opt_distance) from reqm and reqm_hbond
		// reqm: equilibrium internuclear separation 
		//       (sum of the vdW radii of two like atoms (A)) in the case of vdW
		// reqm_hbond: equilibrium internuclear separation
		// 	 (sum of the vdW radii of two like atoms (A)) in the case of hbond 
		float opt_distance;

		uint atom1_type_vdw_hb = KerConstStatic_atom1_types_reqm [atom1_typeid];
     	uint atom2_type_vdw_hb = KerConstStatic_atom2_types_reqm [atom2_typeid];

		if (KerConstStatic_intraE_contributors_const[3*contributor_counter + 2] == 1)	// H-bond
		{
			opt_distance = KerConstStatic_reqm_hbond [atom1_type_vdw_hb] + KerConstStatic_reqm_hbond [atom2_type_vdw_hb];
		}
		else	// Van der Waals
		{
			opt_distance = 0.5f*(KerConstStatic_reqm [atom1_type_vdw_hb] + KerConstStatic_reqm [atom2_type_vdw_hb]);
		}

		// Getting smoothed distance
		// smoothed_distance = function(atomic_distance, opt_distance)
		float smoothed_distance;
		float delta_distance = 0.5f*DockConst_smooth;

		/*printf("delta_distance: %f\n", delta_distance);*/

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

		float smoothed_distance_pow_2 = smoothed_distance*smoothed_distance; 
		float inverse_smoothed_distance_pow_2  = (1.0f / smoothed_distance_pow_2);
		float inverse_smoothed_distance_pow_4  = inverse_smoothed_distance_pow_2 * inverse_smoothed_distance_pow_2;
		float inverse_smoothed_distance_pow_6  = inverse_smoothed_distance_pow_4 * inverse_smoothed_distance_pow_2;
		float inverse_smoothed_distance_pow_10 = inverse_smoothed_distance_pow_6 * inverse_smoothed_distance_pow_4;
		float inverse_smoothed_distance_pow_12 = inverse_smoothed_distance_pow_6 * inverse_smoothed_distance_pow_6;

		// Calculating energy contributions
		// Cuttoff1: internuclear-distance at 8A only for vdw and hbond.
		if (atomic_distance < 8.0f) 
		{
			// Calculating van der Waals / hydrogen bond term
			partialE1 += KerConstStatic_VWpars_AC_const [atom1_typeid*DockConst_num_of_atypes+atom2_typeid]*inverse_smoothed_distance_pow_12;

			float tmp_pE2 = KerConstStatic_VWpars_BD_const [atom1_typeid*DockConst_num_of_atypes+atom2_typeid];

			if (KerConstStatic_intraE_contributors_const[3*contributor_counter + 2] == 1)	// H-bond
				partialE2 -= tmp_pE2 * inverse_smoothed_distance_pow_10;
			else	// Van der Waals
				partialE2 -= tmp_pE2 * inverse_smoothed_distance_pow_6;
		} // if cuttoff1 - internuclear-distance at 8A

		// Calculating energy contributions
		// Cuttoff2: internuclear-distance at 20.48A only for el and sol.
		if (atomic_distance < 20.48f)
		{
			// Calculating electrostatic term
			partialE3 += (  
							(DockConst_coeff_elec*KerConstStatic_atom_charges_const[atom1_id]*KerConstStatic_atom_charges_const[atom2_id]) 
							/
							(
								atomic_distance*(
									-8.5525f + (86.9525f / (1.0f + 7.7839f*exp(-0.3154f*atomic_distance)))
								)
							)       
						);

			// Calculating desolvation term
			partialE4 += (
				  (KerConstStatic_dspars_S_const [atom1_typeid] + DockConst_qasp*fabs(KerConstStatic_atom_charges_const[atom1_id])) * KerConstStatic_dspars_V_const [atom2_typeid] + 
				  (KerConstStatic_dspars_S_const [atom2_typeid] + DockConst_qasp*fabs(KerConstStatic_atom_charges_const[atom2_id])) * KerConstStatic_dspars_V_const [atom1_typeid]) * 
				 DockConst_coeff_desolv*exp(-0.0386f*distance_pow_2);
		} // if cuttoff2 - internuclear-distance at 20.48A
	
		intraE += partialE1 + partialE2 + partialE3 + partialE4;
	
	} // End of contributor_counter for-loop

	// --------------------------------------------------------------
	// Send intramolecular energy to channel
	// --------------------------------------------------------------
	*final_intraE = intraE;

	// --------------------------------------------------------------

	#if defined (DEBUG_KRNL_INTRAE)
	printf("AFTER Out INTRAE CHANNEL\n");
	#endif

	#if defined (DEBUG_ACTIVE_KERNEL)
	printf("	%-20s: %s\n", "Krnl_IntraE", "disabled");
	#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
