#include "auxiliary.h"

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
)
{
#if defined (PRINT_ALL_PC) 
	printf("\n");
	printf("Starting <pose calculation> ... \n");
	printf("\n");
	printf("\t%-40s %u\n", "DockConst_rotbondlist_length: ",	DockConst_rotbondlist_length);
	printf("\t%-40s %u\n", "DockConst_num_of_genes: ",        	DockConst_num_of_genes);
	printf("\t%-40s %u\n", "Host_RunId: ",        			Host_RunId);
#endif

#if defined (ENABLE_TRACE)
	ftrace_region_begin("PC_GENOTYPES_LOOP");
#endif

	float local_genotype[MAX_POPSIZE][ACTUAL_GENOTYPE_LENGTH];

	for (int i = 3; i < DockConst_num_of_genes; i++) {
		for (int j = 0; j < DockConst_pop_size; j++) {
			local_genotype[j][i] = genotype[i][j] * DEG_TO_RAD;
			//local_genotype[i][j] = genotype[i][j] * 0.01745329252;
		}
	}

#if defined (ENABLE_TRACE)
	ftrace_region_end("PC_GENOTYPES_LOOP");
#endif

	float genrot_unitvec[MAX_POPSIZE][3];
	for (int j = 0; j < DockConst_pop_size; j++) {
		float phi = local_genotype[j][3];
		float theta = local_genotype[j][4];

		float sin_theta = sinf(theta);
		float cos_theta = cosf(theta);
		genrot_unitvec[j][0] = sin_theta*cosf(phi);
		genrot_unitvec[j][1] = sin_theta*sinf(phi);
		genrot_unitvec[j][2] = cos_theta;
	}

#if defined (ENABLE_TRACE)
	ftrace_region_begin("PC_MAIN_LOOP");
#endif

	for (int rotation_counter = 0; rotation_counter < DockConst_rotbondlist_length; rotation_counter++) {
		int rotation_list_element = PC_rotlist[rotation_counter];
		if ((rotation_list_element & RLIST_DUMMY_MASK) == 0) {	// If not dummy rotation
			uint atom_id = rotation_list_element & RLIST_ATOMID_MASK;

			for (int j = 0; j < DockConst_pop_size; j++) {

				// Capturing atom coordinates
				float atom_to_rotate[3];

				if ((rotation_list_element & RLIST_FIRSTROT_MASK) != 0)	{ // If first rotation of this atom
					atom_to_rotate[0] = PC_ref_coords_x[atom_id];
					atom_to_rotate[1] = PC_ref_coords_y[atom_id];
					atom_to_rotate[2] = PC_ref_coords_z[atom_id];

				} else {
					atom_to_rotate[0] = local_coords_x[atom_id][j];
					atom_to_rotate[1] = local_coords_y[atom_id][j];
					atom_to_rotate[2] = local_coords_z[atom_id][j];
				}

				// Capturing rotation vectors and angle
				float rotation_unitvec[3];
				float rotation_movingvec[3];
				float rotation_angle;

				if ((rotation_list_element & RLIST_GENROT_MASK) != 0) {	// If general rotation
					rotation_unitvec[0] = genrot_unitvec[j][0];
					rotation_unitvec[1] = genrot_unitvec[j][1];
					rotation_unitvec[2] = genrot_unitvec[j][2];
                  
					rotation_movingvec[0] = genotype[0][j];
					rotation_movingvec[1] = genotype[1][j];
					rotation_movingvec[2] = genotype[2][j];

					rotation_angle = local_genotype[j][5];
				}
				else	// If rotating around rotatable bond
				{
					uint rotbond_id = (rotation_list_element & RLIST_RBONDID_MASK) >> RLIST_RBONDID_SHIFT;

					rotation_unitvec[0] = PC_rotbonds_unit_vectors[3*rotbond_id];
					rotation_unitvec[1] = PC_rotbonds_unit_vectors[3*rotbond_id+1];
					rotation_unitvec[2] = PC_rotbonds_unit_vectors[3*rotbond_id+2];
				
					rotation_movingvec[0] = PC_rotbonds_moving_vectors[3*rotbond_id];
					rotation_movingvec[1] = PC_rotbonds_moving_vectors[3*rotbond_id+1];
					rotation_movingvec[2] = PC_rotbonds_moving_vectors[3*rotbond_id+2];

					rotation_angle = local_genotype[j][6+rotbond_id];

					// In addition performing the first movement 
					// which is needed only if rotating around rotatable bond
					atom_to_rotate[0] -= rotation_movingvec[0];
					atom_to_rotate[1] -= rotation_movingvec[1];
					atom_to_rotate[2] -= rotation_movingvec[2];
				}

				// Performing rotation
				float quatrot_left_q, quatrot_left_x, quatrot_left_y, quatrot_left_z;
				float quatrot_temp_q, quatrot_temp_x, quatrot_temp_y, quatrot_temp_z;

				rotation_angle = rotation_angle * 0.5f;

				float sin_angle, cos_angle;
				sin_angle      = sinf(rotation_angle);
				cos_angle      = cosf(rotation_angle);
				quatrot_left_q = cos_angle;
				quatrot_left_x = sin_angle * rotation_unitvec[0];
				quatrot_left_y = sin_angle * rotation_unitvec[1];
				quatrot_left_z = sin_angle * rotation_unitvec[2];
			

				if ((rotation_list_element & RLIST_GENROT_MASK) != 0)	// If general rotation, 
					// two rotations should be performed 
					// (multiplying the quaternions)
				{
					const float  ref_ori_quats_const_q = PC_ref_orientation_quats[4*Host_RunId];
					const float  ref_ori_quats_const_x = PC_ref_orientation_quats[4*Host_RunId+1];
					const float  ref_ori_quats_const_y = PC_ref_orientation_quats[4*Host_RunId+2];
					const float  ref_ori_quats_const_z = PC_ref_orientation_quats[4*Host_RunId+3];

					// Calculating quatrot_left*ref_orientation_quats_const, 
					// which means that reference orientation rotation is the first
					quatrot_temp_q = quatrot_left_q;
					quatrot_temp_x = quatrot_left_x;
					quatrot_temp_y = quatrot_left_y;
					quatrot_temp_z = quatrot_left_z;

					// Taking the first element of ref_orientation_quats_const member
					quatrot_left_q = esa_dot4_e(       quatrot_temp_q,         quatrot_temp_x,         quatrot_temp_y,         quatrot_temp_z,
								    ref_ori_quats_const_q, -ref_ori_quats_const_x, -ref_ori_quats_const_y, -ref_ori_quats_const_z);
					quatrot_left_x = esa_dot4_e(       quatrot_temp_q,         quatrot_temp_x,         quatrot_temp_y,         quatrot_temp_z,
								    ref_ori_quats_const_x,  ref_ori_quats_const_q,  ref_ori_quats_const_z, -ref_ori_quats_const_y);
					quatrot_left_y = esa_dot4_e(       quatrot_temp_q,         quatrot_temp_x,         quatrot_temp_y,         quatrot_temp_z,
								    ref_ori_quats_const_y, -ref_ori_quats_const_z,  ref_ori_quats_const_q,  ref_ori_quats_const_x);
					quatrot_left_z = esa_dot4_e(       quatrot_temp_q,         quatrot_temp_x,         quatrot_temp_y,         quatrot_temp_z,
								    ref_ori_quats_const_z,  ref_ori_quats_const_y, -ref_ori_quats_const_x,  ref_ori_quats_const_q);
				}

				quatrot_temp_q = esa_dot3_e(  -quatrot_left_x,    -quatrot_left_y,    -quatrot_left_z,
                                                            atom_to_rotate[0],  atom_to_rotate[1],  atom_to_rotate[2]);
				quatrot_temp_x = esa_dot3_e(   quatrot_left_q,    -quatrot_left_z,     quatrot_left_y,
                                                            atom_to_rotate[0],  atom_to_rotate[1],  atom_to_rotate[2]);
				quatrot_temp_y = esa_dot3_e(   quatrot_left_z,     quatrot_left_q,    -quatrot_left_x,
                                                            atom_to_rotate[0],  atom_to_rotate[1],  atom_to_rotate[2]);
				quatrot_temp_z = esa_dot3_e(  -quatrot_left_y,     quatrot_left_x,     quatrot_left_q,
                                                            atom_to_rotate[0],  atom_to_rotate[1],  atom_to_rotate[2]);

				atom_to_rotate[0] = esa_dot4_e( quatrot_temp_q, quatrot_temp_x, quatrot_temp_y, quatrot_temp_z,
                                                               -quatrot_left_x, quatrot_left_q,-quatrot_left_z, quatrot_left_y);
				atom_to_rotate[1] = esa_dot4_e( quatrot_temp_q, quatrot_temp_x, quatrot_temp_y, quatrot_temp_z,
                                                               -quatrot_left_y, quatrot_left_z, quatrot_left_q,-quatrot_left_x);
				atom_to_rotate[2] = esa_dot4_e( quatrot_temp_q, quatrot_temp_x, quatrot_temp_y, quatrot_temp_z,
							       -quatrot_left_z,-quatrot_left_y, quatrot_left_x, quatrot_left_q);
				// Performing final movement and storing values
				local_coords_x[atom_id][j] = (float)(atom_to_rotate[0] + rotation_movingvec[0]);
				local_coords_y[atom_id][j] = (float)(atom_to_rotate[1] + rotation_movingvec[1]);
				local_coords_z[atom_id][j] = (float)(atom_to_rotate[2] + rotation_movingvec[2]);

			} // Loop over j (individuals)

		} // End if-statement not dummy rotation

	} // End rotation_counter for-loop
    
#if defined (ENABLE_TRACE)
	ftrace_region_end("PC_MAIN_LOOP");
#endif

#if defined (PRINT_ALL_PC) 
	printf("\n");
	printf("Finishing <pose calculation>\n");
	printf("\n");
#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
