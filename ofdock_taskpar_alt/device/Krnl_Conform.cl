#include "auxiliary_conform.cl"

// --------------------------------------------------------------------------
// The function changes the conformation of myligand according to 
// the genotype given by the second parameter.
// Originally from: processligand.c
// --------------------------------------------------------------------------
__kernel __attribute__ ((reqd_work_group_size(1,1,1)))
/*
void Krnl_Conform(__global const int*             restrict GlobLigand_rotbonds,
		  __global const char*             restrict GlobLigand_atom_rotbonds,
		  __global const float*           restrict GlobLigand_rotbonds_moving_vectors,
		  __global const float*           restrict GlobLigand_rotbonds_unit_vectors,
		  __global const Ligandconstant*  restrict LigConst)
*/
void Krnl_Conform(
             __global const float*           restrict GlobFgrids,
	     __global       float*           restrict GlobPopulationCurrent,
	     __global       float*           restrict GlobEnergyCurrent,
	     __global       float*           restrict GlobPopulationNext,
	     __global       float*           restrict GlobEnergyNext,
             __global const float*           restrict GlobPRNG,
	     __global const kernelconstant* restrict KerConst,
	     __global const Dockparameters* restrict DockConst)
{

#ifdef EMULATOR
	printf("Krnl Conform!!\n");
#endif
	// --------------------------------------------------------------
	// Wait for ligand data and genotypes
	// --------------------------------------------------------------
/*
	__local float myligand_atom_idxyzq[MAX_NUM_OF_ATOMS*5];
*/

	__local float loc_coords_x[MAX_NUM_OF_ATOMS];
	__local float loc_coords_y[MAX_NUM_OF_ATOMS];
	__local float loc_coords_z[MAX_NUM_OF_ATOMS];
	__local float genotype[ACTUAL_GENOTYPE_LENGTH];

	for (uint pipe_cnt=0; pipe_cnt<DockConst->num_of_atoms; pipe_cnt++) {
		loc_coords_x[pipe_cnt] = read_channel_altera(chan_GA2Conf_x);
		loc_coords_y[pipe_cnt] = read_channel_altera(chan_GA2Conf_y);
		loc_coords_z[pipe_cnt] = read_channel_altera(chan_GA2Conf_z);
	}


	for (uint pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
		genotype[pipe_cnt] = read_channel_altera(chan_GA2Conf_genotype);
	}
	// --------------------------------------------------------------






	const float genrot_movvec [3] = {0.0f,0.0f,0.0f};
	float genrot_unitvec  [3];
	float movvec_to_origo [3];

	float phi, theta;
	uint atom_id, rotbond_id;

	phi   = genotype [3]*DEG_TO_RAD;
	theta = genotype [4]*DEG_TO_RAD;

	#if defined (NATIVE_PRECISION)
	genrot_unitvec [0] = native_sin(theta)*native_cos(phi);
	genrot_unitvec [1] = native_sin(theta)*native_sin(phi);
	genrot_unitvec [2] = native_cos(theta);
	#else
	genrot_unitvec [0] = sin(theta)*cos(phi);
	genrot_unitvec [1] = sin(theta)*sin(phi);
	genrot_unitvec [2] = cos(theta);
	#endif
	
// NOT SURE IF THIS MUST BE KEPT WITH THE 
// SINGLE FOR-LOOP/DATA PARALLEL IMPLEMENTATION
/*
	//moving ligand to origo
	get_movvec_to_origo(DockConst->num_of_atoms, loc_coords_x, loc_coords_y, loc_coords_z, movvec_to_origo);
	move_ligand        (DockConst->num_of_atoms, loc_coords_x, loc_coords_y, loc_coords_z, movvec_to_origo);
*/
	int rotation_list_element;
	float atom_to_rotate[3];
	float rotation_angle;
	float genrotangle;
	float sin_angle;
	float rotation_unitvec[3], rotation_movingvec[3];
	float quatrot_left_x, quatrot_left_y, quatrot_left_z, quatrot_left_q;
	float quatrot_temp_x, quatrot_temp_y, quatrot_temp_z, quatrot_temp_q;

	genrotangle = genotype[5]*DEG_TO_RAD;

	// **********************************************
	// ADD VENDOR SPECIFIC PRAGMA
	// **********************************************
	LOOP_CHANGE_CONFORM_1:
	for (uint rotation_counter = 0; rotation_counter < DockConst->rotbondlist_length; rotation_counter++)
	{
		rotation_list_element = KerConst->rotlist_const[rotation_counter];

		if ((rotation_list_element & RLIST_DUMMY_MASK) == 0)	//if not dummy rotation
		{
			atom_id = rotation_list_element & RLIST_ATOMID_MASK;

			//capturing atom coordinates
			if ((rotation_list_element & RLIST_FIRSTROT_MASK) != 0)	//if first rotation of this atom
			{
				atom_to_rotate[0] = KerConst->ref_coords_x_const[atom_id];
				atom_to_rotate[1] = KerConst->ref_coords_y_const[atom_id];
				atom_to_rotate[2] = KerConst->ref_coords_z_const[atom_id];
			}
			else
			{
				atom_to_rotate[0] = loc_coords_x[atom_id];
				atom_to_rotate[1] = loc_coords_y[atom_id];
				atom_to_rotate[2] = loc_coords_z[atom_id];
			}

			//capturing rotation vectors and angle
			if ((rotation_list_element & RLIST_GENROT_MASK) != 0)	//if general rotation
			{
				rotation_unitvec[0] = genrot_unitvec[0];
				rotation_unitvec[1] = genrot_unitvec[1];
				rotation_unitvec[2] = genrot_unitvec[2];

				rotation_angle = genrotangle;

				rotation_movingvec[0] = genotype[0];
				rotation_movingvec[1] = genotype[1];
				rotation_movingvec[2] = genotype[2];
			}
			else	//if rotating around rotatable bond
			{
				rotbond_id = (rotation_list_element & RLIST_RBONDID_MASK) >> RLIST_RBONDID_SHIFT;
	
				rotation_unitvec[0] = KerConst->rotbonds_unit_vectors_const[3*rotbond_id];
				rotation_unitvec[1] = KerConst->rotbonds_unit_vectors_const[3*rotbond_id+1];
				rotation_unitvec[2] = KerConst->rotbonds_unit_vectors_const[3*rotbond_id+2];
				rotation_angle = genotype[6+rotbond_id]*DEG_TO_RAD;

				rotation_movingvec[0] = KerConst->rotbonds_moving_vectors_const[3*rotbond_id];
				rotation_movingvec[1] = KerConst->rotbonds_moving_vectors_const[3*rotbond_id+1];
				rotation_movingvec[2] = KerConst->rotbonds_moving_vectors_const[3*rotbond_id+2];

				//in addition performing the first movement 
				//which is needed only if rotating around rotatable bond
				atom_to_rotate[0] -= rotation_movingvec[0];
				atom_to_rotate[1] -= rotation_movingvec[1];
				atom_to_rotate[2] -= rotation_movingvec[2];
			}

			//performing rotation
			rotation_angle = rotation_angle/2;
			quatrot_left_q = cos(rotation_angle);
			sin_angle = sin(rotation_angle);

			quatrot_left_x = sin_angle*rotation_unitvec[0];
			quatrot_left_y = sin_angle*rotation_unitvec[1];
			quatrot_left_z = sin_angle*rotation_unitvec[2];

			if ((rotation_list_element & RLIST_GENROT_MASK) != 0)	//if general rotation, 
										//two rotations should be performed 
										//(multiplying the quaternions)
			{
				//calculating quatrot_left*ref_orientation_quats_const, 
				//which means that reference orientation rotation is the first
				quatrot_temp_q = quatrot_left_q;
				quatrot_temp_x = quatrot_left_x;
				quatrot_temp_y = quatrot_left_y;
				quatrot_temp_z = quatrot_left_z;

				// L30nardoSV: taking the first one
				quatrot_left_q = quatrot_temp_q*KerConst->ref_orientation_quats_const[0]-
						 quatrot_temp_x*KerConst->ref_orientation_quats_const[1]-
						 quatrot_temp_y*KerConst->ref_orientation_quats_const[2]-
						 quatrot_temp_z*KerConst->ref_orientation_quats_const[3];
				quatrot_left_x = quatrot_temp_q*KerConst->ref_orientation_quats_const[1]+
						 KerConst->ref_orientation_quats_const[0]*quatrot_temp_x+
						 quatrot_temp_y*KerConst->ref_orientation_quats_const[3]-
						 KerConst->ref_orientation_quats_const[2]*quatrot_temp_z;
				quatrot_left_y = quatrot_temp_q*KerConst->ref_orientation_quats_const[2]+
						 KerConst->ref_orientation_quats_const[0]*quatrot_temp_y+
						 KerConst->ref_orientation_quats_const[1]*quatrot_temp_z-
						 quatrot_temp_x*KerConst->ref_orientation_quats_const[3];
				quatrot_left_z = quatrot_temp_q*KerConst->ref_orientation_quats_const[3]+
						 KerConst->ref_orientation_quats_const[0]*quatrot_temp_z+
						 quatrot_temp_x*KerConst->ref_orientation_quats_const[2]-
						 KerConst->ref_orientation_quats_const[1]*quatrot_temp_y;

			}

			quatrot_temp_q = 0 -
					 quatrot_left_x*atom_to_rotate [0] -
					 quatrot_left_y*atom_to_rotate [1] -
					 quatrot_left_z*atom_to_rotate [2];
			quatrot_temp_x = quatrot_left_q*atom_to_rotate [0] +
					 quatrot_left_y*atom_to_rotate [2] -
					 quatrot_left_z*atom_to_rotate [1];
			quatrot_temp_y = quatrot_left_q*atom_to_rotate [1] -
					 quatrot_left_x*atom_to_rotate [2] +
					 quatrot_left_z*atom_to_rotate [0];
			quatrot_temp_z = quatrot_left_q*atom_to_rotate [2] +
					 quatrot_left_x*atom_to_rotate [1] -
					 quatrot_left_y*atom_to_rotate [0];

			atom_to_rotate [0] = 0 -
					     quatrot_temp_q*quatrot_left_x +
					     quatrot_temp_x*quatrot_left_q -
					     quatrot_temp_y*quatrot_left_z +
					     quatrot_temp_z*quatrot_left_y;
			atom_to_rotate [1] = 0 -
					     quatrot_temp_q*quatrot_left_y +
					     quatrot_temp_x*quatrot_left_z +
					     quatrot_temp_y*quatrot_left_q -
					     quatrot_temp_z*quatrot_left_x;
			atom_to_rotate [2] = 0 -
					     quatrot_temp_q*quatrot_left_z -
					     quatrot_temp_x*quatrot_left_y +
					     quatrot_temp_y*quatrot_left_x +
					     quatrot_temp_z*quatrot_left_q;

			//performing final movement and storing values
			loc_coords_x[atom_id] = atom_to_rotate [0] + rotation_movingvec[0];
			loc_coords_y[atom_id] = atom_to_rotate [1] + rotation_movingvec[1];
			loc_coords_z[atom_id] = atom_to_rotate [2] + rotation_movingvec[2];

		} // End if-statement not dummy rotation
	} // End rotation_counter for-loop



// NOT SURE IF THIS MUST BE KEPT WITH THE 
// SINGLE FOR-LOOP/DATA PARALLEL IMPLEMENTATION
/*
	float genotype_copy [3];
	genotype_copy [0] = genotype [0];
	genotype_copy [1] = genotype [1];
	genotype_copy [2] = genotype [2];
	
	move_ligand(LigConst->num_of_atoms, myligand_atom_idxyzq, genotype_copy);
	
	#if defined (DEBUG_KERNEL_CHANGE_CONFORM)
	for (atom_id=0; atom_id<LigConst->num_of_atoms; atom_id++) {
		printf("Moved point (final values) (x,y,z): %f, %f, %f\n", 
		       myligand_atom_idxyzq [atom_id*5+1],
		       myligand_atom_idxyzq [atom_id*5+2],
		       myligand_atom_idxyzq [atom_id*5+3]);
	}
	#endif

*/


	// --------------------------------------------------------------
	// Send ligand data to InterE and IntraE Kernels
	// --------------------------------------------------------------
	for (uint pipe_cnt=0; pipe_cnt<DockConst->num_of_atoms; pipe_cnt++) {
		write_channel_altera(chan_Conf2Intere_x, loc_coords_x[pipe_cnt]);
		write_channel_altera(chan_Conf2Intere_y, loc_coords_y[pipe_cnt]);
		write_channel_altera(chan_Conf2Intere_z, loc_coords_z[pipe_cnt]);

		write_channel_altera(chan_Conf2Intrae_x, loc_coords_x[pipe_cnt]);
		write_channel_altera(chan_Conf2Intrae_y, loc_coords_y[pipe_cnt]);
		write_channel_altera(chan_Conf2Intrae_z, loc_coords_z[pipe_cnt]);
	}
	// --------------------------------------------------------------


}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
