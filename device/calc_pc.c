#include "auxiliary.c"
#include "math.h"

// --------------------------------------------------------------------------
// Conform changes the conformation of the ligand according to 
// the genotype fed by any producer logic/kernel (IC, GG, LSs).
// Originally from: processligand.c
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
			ushort				Host_RunId,

	const	float*	restrict	genotype,
			float*	restrict 	local_coords_x,
			float*	restrict	local_coords_y,
			float*	restrict	local_coords_z
)
{
	#if defined (DEBUG_KRNL_Conform) 
	printf("\n");
	printf("%-40s %u\n", "DockConst->rotbondlist_length: ",	DockConst_rotbondlist_length);
	printf("%-40s %u\n", "DockConst_num_of_genes: ",        DockConst_num_of_genes);	
	#endif

	// Check best practices guide
	// Table 11. Effects of numbanks and bankwidth on the Bank Geometry ...
	// Only first three indexes of the lower array are used
	// however size of lower array was declared as 4, 
	// just to keep sizes equal to power of 2
	// __local float  __attribute__((numbanks(8), bankwidth(16))) local_coords[MAX_NUM_OF_ATOMS][4];

	int rotlist_localcache [MAX_NUM_OF_ROTATIONS];
	for (ushort c = 0; c < DockConst_rotbondlist_length; c++) {
		rotlist_localcache [c] = PC_rotlist [c];
	}

	float local_genotype [ACTUAL_GENOTYPE_LENGTH];
	for (uchar i=0; i<DockConst_num_of_genes; i++) {
		if (i < 3) {
			local_genotype [i] = genotype[i];
		} else {
			local_genotype [i] = genotype[i] * DEG_TO_RAD;
		}
	}

	float phi = local_genotype[3];
	float theta = local_genotype[4];
	float genrotangle = local_genotype[5];

	#if defined (DEBUG_ACTIVE_KERNEL)
	if (active == 0x00) {printf("	%-20s: %s\n", "Krnl_Conform", "must be disabled");}
	#endif

	for (ushort rotation_counter = 0; rotation_counter < DockConst_rotbondlist_length; rotation_counter++)
	{
		int rotation_list_element = rotlist_localcache [rotation_counter];

		if ((rotation_list_element & RLIST_DUMMY_MASK) == 0)	// If not dummy rotation
		{
			uint atom_id = rotation_list_element & RLIST_ATOMID_MASK;

			// Capturing atom coordinates
			float atom_to_rotate[3];

			if ((rotation_list_element & RLIST_FIRSTROT_MASK) != 0)	// If first rotation of this atom
			{	
				atom_to_rotate[0] = PC_ref_coords_x [atom_id];
				atom_to_rotate[1] = PC_ref_coords_y [atom_id];
				atom_to_rotate[2] = PC_ref_coords_z [atom_id];
			}
			else
			{	
				atom_to_rotate[0] = local_coords_x[atom_id];
				atom_to_rotate[1] = local_coords_y[atom_id];
				atom_to_rotate[2] = local_coords_z[atom_id];
			}

			// Capturing rotation vectors and angle
			float rotation_unitvec[3];
			float rotation_movingvec[3];
			float rotation_angle;

			if ((rotation_list_element & RLIST_GENROT_MASK) != 0)	// If general rotation
			{
				float  sin_theta, cos_theta;
				float genrot_unitvec[3];
				sin_theta = sin(theta);
				cos_theta = cos(theta);
				genrot_unitvec[0] = sin_theta*cos(phi);
				genrot_unitvec[1] = sin_theta*sin(phi);
				genrot_unitvec[2] = cos_theta;

				rotation_unitvec[0] = genrot_unitvec[0];
				rotation_unitvec[1] = genrot_unitvec[1];
				rotation_unitvec[2] = genrot_unitvec[2];
				
				rotation_movingvec[0] = local_genotype[0];
				rotation_movingvec[1] = local_genotype[1];
				rotation_movingvec[2] = local_genotype[2];

				rotation_angle = genrotangle;
			}
			else	// If rotating around rotatable bond
			{
				uint rotbond_id = (rotation_list_element & RLIST_RBONDID_MASK) >> RLIST_RBONDID_SHIFT;

				rotation_unitvec[0] = PC_rotbonds_unit_vectors[rotbond_id];
				rotation_unitvec[1] = PC_rotbonds_unit_vectors[rotbond_id];
				rotation_unitvec[2] = PC_rotbonds_unit_vectors[rotbond_id];
				
				rotation_movingvec[0] = PC_rotbonds_moving_vectors[rotbond_id];
				rotation_movingvec[1] = PC_rotbonds_moving_vectors[rotbond_id];
				rotation_movingvec[2] = PC_rotbonds_moving_vectors[rotbond_id];

				rotation_angle = local_genotype[6+rotbond_id];

				// In addition performing the first movement 
				// which is needed only if rotating around rotatable bond
				atom_to_rotate[0] -= rotation_movingvec[0];
				atom_to_rotate[1] -= rotation_movingvec[1];
				atom_to_rotate[2] -= rotation_movingvec[2];
			}

			// Performing rotation
			float quatrot_left_q, quatrot_left_x, quatrot_left_y, quatrot_left_z;
			float quatrot_temp_q, quatrot_temp_x, quatrot_temp_y, quatrot_temp_z;

			rotation_angle = rotation_angle*0.5f;

			float sin_angle, cos_angle;
			sin_angle      = sin(rotation_angle);
			cos_angle      = cos(rotation_angle);
			quatrot_left_q = cos_angle;
			quatrot_left_x = sin_angle * rotation_unitvec[0];
			quatrot_left_y = sin_angle * rotation_unitvec[1];
			quatrot_left_z = sin_angle * rotation_unitvec[2];
			

			if ((rotation_list_element & RLIST_GENROT_MASK) != 0)	// If general rotation, 
										// two rotations should be performed 
										// (multiplying the quaternions)
			{
				const float  ref_ori_quats_const_q = PC_ref_orientation_quats[Host_RunId];
				const float  ref_ori_quats_const_x = PC_ref_orientation_quats[Host_RunId];
				const float  ref_ori_quats_const_y = PC_ref_orientation_quats[Host_RunId];
				const float  ref_ori_quats_const_z = PC_ref_orientation_quats[Host_RunId];

				// Calculating quatrot_left*ref_orientation_quats_const, 
				// which means that reference orientation rotation is the first
				quatrot_temp_q = quatrot_left_q;
				quatrot_temp_x = quatrot_left_x;
				quatrot_temp_y = quatrot_left_y;
				quatrot_temp_z = quatrot_left_z;
				float quatrot_temp[4] = { quatrot_temp_q, quatrot_temp_x, quatrot_temp_y, quatrot_temp_z };

				// Taking the first element of ref_orientation_quats_const member
				float ref_4q[4] = {  ref_ori_quats_const_q, -ref_ori_quats_const_x, -ref_ori_quats_const_y, -ref_ori_quats_const_z };
				float ref_4x[4] = {  ref_ori_quats_const_x,  ref_ori_quats_const_q, -ref_ori_quats_const_z,  ref_ori_quats_const_y };
				float ref_4y[4] = { -ref_ori_quats_const_y,  ref_ori_quats_const_z,  ref_ori_quats_const_q,  ref_ori_quats_const_x };
				float ref_4z[4] = {  ref_ori_quats_const_z, -ref_ori_quats_const_y,  ref_ori_quats_const_x,  ref_ori_quats_const_q };

				quatrot_left_q = esa_dot4(quatrot_temp, ref_4q);
				quatrot_left_x = esa_dot4(quatrot_temp, ref_4x);
				quatrot_left_y = esa_dot4(quatrot_temp, ref_4y);
				quatrot_left_z = esa_dot4(quatrot_temp, ref_4z);
			}

			float left_3q[3] = {-quatrot_left_x, -quatrot_left_y, -quatrot_left_z };
			float left_3x[3] = { quatrot_left_q, -quatrot_left_z,  quatrot_left_y };
			float left_3y[3] = { quatrot_left_z,  quatrot_left_q, -quatrot_left_x };
			float left_3z[3] = {-quatrot_left_y,  quatrot_left_x,  quatrot_left_q };
			
			quatrot_temp_q = esa_dot3(left_3q, atom_to_rotate);
			quatrot_temp_x = esa_dot3(left_3x, atom_to_rotate);
			quatrot_temp_y = esa_dot3(left_3y, atom_to_rotate);
			quatrot_temp_z = esa_dot3(left_3z, atom_to_rotate);
			float quatrot_temp_2[4] = { quatrot_temp_q, quatrot_temp_x, quatrot_temp_y, quatrot_temp_z };

			float left_4x[4] = { -quatrot_left_x,  quatrot_left_q, -quatrot_left_z,  quatrot_left_y };
			float left_4y[4] = { -quatrot_left_y,  quatrot_left_z,  quatrot_left_q, -quatrot_left_x };
			float left_4z[4] = { -quatrot_left_z, -quatrot_left_y,  quatrot_left_x,  quatrot_left_q };

			atom_to_rotate[0] = esa_dot4(quatrot_temp_2, left_4x);
			atom_to_rotate[1] = esa_dot4(quatrot_temp_2, left_4y);
			atom_to_rotate[2] = esa_dot4(quatrot_temp_2, left_4z);

			// Performing final movement and storing values
			local_coords_x[atom_id] = atom_to_rotate[0] + rotation_movingvec[0];
			local_coords_y[atom_id] = atom_to_rotate[1] + rotation_movingvec[1];
			local_coords_z[atom_id] = atom_to_rotate[2] + rotation_movingvec[2];
		} // End if-statement not dummy rotation

		//mem_fence(CLK_LOCAL_MEM_FENCE);

	} // End rotation_counter for-loop

#if defined (DEBUG_ACTIVE_KERNEL)
printf("	%-20s: %s\n", "Krnl_Conform", "disabled");
#endif

}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
