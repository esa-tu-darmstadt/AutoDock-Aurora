#include "calcenergy.h"

int prepare_conststatic_fields_for_aurora(Liganddata* 	       myligand_reference,
				 	Dockpars*   	       mypars,
				 	float*      	       cpu_ref_ori_angles,
				 	kernelconstant_static* KerConstStatic)
{
	int i, j;
	int type_id1, type_id2;
	float* floatpoi;
	char* charpoi;
	//float phi, theta, genrotangle;

	// arrays to store intermmediate data
	float atom_charges[MAX_NUM_OF_ATOMS];
	char  atom_types[MAX_NUM_OF_ATOMS];
	char  intraE_contributors[3*MAX_INTRAE_CONTRIBUTORS];

	float reqm [ATYPE_NUM];
	float reqm_hbond [ATYPE_NUM];
	unsigned int atom1_types_reqm [ATYPE_NUM];
	unsigned int atom2_types_reqm [ATYPE_NUM];

	float VWpars_AC[MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
	float VWpars_BD[MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
	float dspars_S[MAX_NUM_OF_ATYPES];
	float dspars_V[MAX_NUM_OF_ATYPES];
	int rotlist[MAX_NUM_OF_ROTATIONS];

	//charges and type id-s
	floatpoi = atom_charges;
	charpoi = atom_types;

	for (i=0; i < myligand_reference->num_of_atoms; i++)
	{
		*floatpoi = (float) myligand_reference->atom_idxyzq[i][4];
		*charpoi = (char) myligand_reference->atom_idxyzq[i][0];
		floatpoi++;
		charpoi++;
	}

	//intramolecular energy contributors
	myligand_reference->num_of_intraE_contributors = 0;
	for (i=0; i<myligand_reference->num_of_atoms-1; i++)
		for (j=i+1; j<myligand_reference->num_of_atoms; j++)
		{
			if (myligand_reference->intraE_contributors[i][j])
				myligand_reference->num_of_intraE_contributors++;
		}

	if (myligand_reference->num_of_intraE_contributors > MAX_INTRAE_CONTRIBUTORS)
	{
		printf("Error: number of intramolecular energy contributor is too high!\n");
		fflush(stdout);
		return 1;
	}

	charpoi = intraE_contributors;
	for (i=0; i<myligand_reference->num_of_atoms-1; i++)
		for (j=i+1; j<myligand_reference->num_of_atoms; j++)
		{
			if (myligand_reference->intraE_contributors[i][j] == 1)
			{
				*charpoi = (char) i;
				charpoi++;
				*charpoi = (char) j;
				charpoi++;

				type_id1 = (int) myligand_reference->atom_idxyzq [i][0];
				type_id2 = (int) myligand_reference->atom_idxyzq [j][0];
				if (is_H_bond(myligand_reference->atom_types[type_id1], myligand_reference->atom_types[type_id2]) != 0)
					*charpoi = (char) 1;
				else
					*charpoi = (char) 0;
				charpoi++;
			}
		}

    // -------------------------------------------
    // Smoothed pairwise potentials
    // -------------------------------------------
	// reqm, reqm_hbond: equilibrium internuclear separation for vdW and hbond
	for (i= 0; i<ATYPE_NUM/*myligand_reference->num_of_atypes*/; i++) {
		reqm[i]       = myligand_reference->reqm[i];
		reqm_hbond[i] = myligand_reference->reqm_hbond[i];

		atom1_types_reqm [i] = myligand_reference->atom1_types_reqm[i];
        	atom2_types_reqm [i] = myligand_reference->atom2_types_reqm[i];
	}
	// -------------------------------------------

	//van der Waals parameters
	for (i=0; i<myligand_reference->num_of_atypes; i++)
		for (j=0; j<myligand_reference->num_of_atypes; j++)
		{
			if (is_H_bond(myligand_reference->atom_types[i], myligand_reference->atom_types[j]) != 0)
			{
				floatpoi = VWpars_AC + i*myligand_reference->num_of_atypes + j;
				*floatpoi = (float) myligand_reference->VWpars_C[i][j];
				floatpoi = VWpars_AC + j*myligand_reference->num_of_atypes + i;
				*floatpoi = (float) myligand_reference->VWpars_C[j][i];

				floatpoi = VWpars_BD + i*myligand_reference->num_of_atypes + j;
				*floatpoi = (float) myligand_reference->VWpars_D[i][j];
				floatpoi = VWpars_BD + j*myligand_reference->num_of_atypes + i;
				*floatpoi = (float) myligand_reference->VWpars_D[j][i];
			}
			else
			{
				floatpoi = VWpars_AC + i*myligand_reference->num_of_atypes + j;
				*floatpoi = (float) myligand_reference->VWpars_A[i][j];
				floatpoi = VWpars_AC + j*myligand_reference->num_of_atypes + i;
				*floatpoi = (float) myligand_reference->VWpars_A[j][i];

				floatpoi = VWpars_BD + i*myligand_reference->num_of_atypes + j;
				*floatpoi = (float) myligand_reference->VWpars_B[i][j];
				floatpoi = VWpars_BD + j*myligand_reference->num_of_atypes + i;
				*floatpoi = (float) myligand_reference->VWpars_B[j][i];
			}
		}

	//desolvation parameters
	for (i=0; i<myligand_reference->num_of_atypes; i++)
	{
		dspars_S[i] = myligand_reference->solpar[i];
		dspars_V[i] = myligand_reference->volume[i];
	}

	//generate rotation list
	if (gen_rotlist(
			myligand_reference, rotlist
                        ) != 0)
	{
		printf("Error: number of required rotations is too high!\n");
		return 1;
	}

	int m;
	for (m=0;m<MAX_NUM_OF_ATOMS;m++) {
		KerConstStatic->atom_charges_const[m] = atom_charges[m]; 
	}

	for (m=0;m<MAX_NUM_OF_ATOMS;m++) { 
		KerConstStatic->atom_types_const[m] = atom_types[m];   
	}

	for (m=0;m<3*MAX_INTRAE_CONTRIBUTORS;m++) { 
		KerConstStatic->intraE_contributors_const[m] = intraE_contributors[m];
	}

	// -------------------------------------------
	// Smoothed pairwise potentials
	// -------------------------------------------
	for (m=0;m<ATYPE_NUM;m++) { KerConstStatic->reqm_const[m]             = reqm[m]; }
	for (m=0;m<ATYPE_NUM;m++) { KerConstStatic->reqm_hbond_const[m]       = reqm_hbond[m]; }
	for (m=0;m<ATYPE_NUM;m++) { KerConstStatic->atom1_types_reqm_const[m] = atom1_types_reqm[m]; }
	for (m=0;m<ATYPE_NUM;m++) { KerConstStatic->atom2_types_reqm_const[m] = atom2_types_reqm[m]; }
	// -------------------------------------------

	for (m=0;m<MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES;m++){ KerConstStatic->VWpars_AC_const[m] = VWpars_AC[m]; }
	for (m=0;m<MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES;m++){ KerConstStatic->VWpars_BD_const[m] = VWpars_BD[m]; }
	for (m=0;m<MAX_NUM_OF_ATYPES;m++) { KerConstStatic->dspars_S_const[m] = dspars_S[m]; }
	for (m=0;m<MAX_NUM_OF_ATYPES;m++) { KerConstStatic->dspars_V_const[m] = dspars_V[m]; }
	for (m=0; m<MAX_NUM_OF_ROTATIONS; m++) {
		KerConstStatic->rotlist_const[m] = rotlist[m];
	}

	//coordinates of reference ligand
	for (i=0; i < myligand_reference->num_of_atoms; i++) {
		KerConstStatic->ref_coords_x_const[i] = myligand_reference->atom_idxyzq[i][1];
		KerConstStatic->ref_coords_y_const[i] = myligand_reference->atom_idxyzq[i][2];
		KerConstStatic->ref_coords_z_const[i] = myligand_reference->atom_idxyzq[i][3];
	}

	//rotatable bond vectors
	for (i=0; i < myligand_reference->num_of_rotbonds; i++) {
		for (j=0; j<3; j++) {
			KerConstStatic->rotbonds_moving_vectors_const[3*i+j] = myligand_reference->rotbonds_moving_vectors[i][j];
			KerConstStatic->rotbonds_unit_vectors_const[3*i+j] = myligand_reference->rotbonds_unit_vectors[i][j];
		}
	}		

	float phi, theta, genrotangle;

	//reference orientation quaternions
	for (unsigned int i=0; i<mypars->num_of_runs; i++)
	{
		//printf("Pregenerated angles for run %d: %f %f %f\n", i, cpu_ref_ori_angles[3*i], cpu_ref_ori_angles[3*i+1], cpu_ref_ori_angles[3*i+2]);
		phi = cpu_ref_ori_angles[3*i]*DEG_TO_RAD;
		theta = cpu_ref_ori_angles[3*i+1]*DEG_TO_RAD;
		genrotangle = cpu_ref_ori_angles[3*i+2]*DEG_TO_RAD;

		KerConstStatic->ref_orientation_quats_const[4*i]   = cosf(genrotangle/2.0f);				//q
		KerConstStatic->ref_orientation_quats_const[4*i+1] = sinf(genrotangle/2.0f)*sinf(theta)*cosf(phi);	//x
		KerConstStatic->ref_orientation_quats_const[4*i+2] = sinf(genrotangle/2.0f)*sinf(theta)*sinf(phi);	//y
		KerConstStatic->ref_orientation_quats_const[4*i+3] = sinf(genrotangle/2.0f)*cosf(theta);		//z
	}

	return 0;
}

void make_reqrot_ordering(char number_of_req_rotations[MAX_NUM_OF_ATOMS],
			  char atom_id_of_numrots[MAX_NUM_OF_ATOMS],
			  int  num_of_atoms)
//The function puts the first array into a descending order and
//performs the same operations on the second array (since element i of
//number_or_req_rotations and element i of atom_id_of_numrots correspond to each other).
//Element i of the former array stores how many rotations have to be perfomed on the atom
//whose atom ID is stored by element i of the latter one. The third parameter has to be equal
//to the number of ligand atoms
{
	int i, j;
	char temp;

	for (j=0; j<num_of_atoms-1; j++)
		for (i=num_of_atoms-2; i>=j; i--)
			if (number_of_req_rotations[i+1] > number_of_req_rotations[i])
			{
				temp = number_of_req_rotations[i];
				number_of_req_rotations[i] = number_of_req_rotations[i+1];
				number_of_req_rotations[i+1] = temp;

				temp = atom_id_of_numrots[i];
				atom_id_of_numrots[i] = atom_id_of_numrots[i+1];
				atom_id_of_numrots[i+1] = temp;
			}

/*	printf("\n\nRotation priority list after re-ordering:\n");
	for (i=0; i<num_of_atoms; i++)
		printf("Rotation of %d (required rots remaining: %d)\n", atom_id_of_numrots[i], number_of_req_rotations[i]);
	printf("\n\n");*/
}

int gen_rotlist(
	Liganddata*		myligand,
	int				rotlist[MAX_NUM_OF_ROTATIONS]
)
// The function generates the rotation list that will be stored in the memory field rotlist_const by
// prepare_conststatic_fields_for_aurora(). The structure of this array is described in that function.
{
	int atom_id, rotb_id, parallel_rot_id, rotlist_id;
	char number_of_req_rotations[MAX_NUM_OF_ATOMS];
	char atom_id_of_numrots[MAX_NUM_OF_ATOMS];
	char atom_wasnt_rotated_yet[MAX_NUM_OF_ATOMS];
	int new_rotlist_element;
	char rotbond_found;
	char rotbond_candidate;
	char remaining_rots_around_rotbonds;

	// Subrotation lists
	// At some point, some code was proposed to build subrotations lists.
	// This was removed thereafter due to lower benefits in performance
	// compared to outer-loop pushing.
	// The printfs below are just vestige from that proposal,
	// which might be retaken eventually.
/*
	char number_of_req_rotations_copy[MAX_NUM_OF_ATOMS];
*/

	myligand->num_of_rotcyc = 0;
	myligand->num_of_rotations_required = 0;

	// On the SX-Aurora machine: kernels are single threaded
	const int num_of_threads_per_block = 1;

	for (atom_id=0; atom_id<num_of_threads_per_block; atom_id++)	//handling special case when num_of_atoms<num_of_threads_per_block
		number_of_req_rotations[atom_id] = 0;

	for (atom_id=0; atom_id<myligand->num_of_atoms; atom_id++)
	{
		atom_id_of_numrots[atom_id] = atom_id;
		atom_wasnt_rotated_yet[atom_id] = 1;

		number_of_req_rotations[atom_id] = 1;

		for (rotb_id=0; rotb_id<myligand->num_of_rotbonds; rotb_id++)
			if (myligand->atom_rotbonds[atom_id][rotb_id] != 0)
				(number_of_req_rotations[atom_id])++;

		myligand->num_of_rotations_required += number_of_req_rotations[atom_id];
	}

	// Subrotation lists
	// At some point, some code was proposed to build subrotations lists.
	// This was removed thereafter due to lower benefits in performance
	// compared to outer-loop pushing.
	// The printfs below are just vestige from that proposal,
	// which might be retaken eventually.
/*
	for (atom_id=0; atom_id<myligand->num_of_atoms; atom_id++)
	{
		number_of_req_rotations_copy[atom_id] = number_of_req_rotations[atom_id];
	}
*/
	rotlist_id=0;
	make_reqrot_ordering(number_of_req_rotations, atom_id_of_numrots, myligand->num_of_atoms);
	while (number_of_req_rotations[0] != 0)	//if the atom with the most remaining rotations has to be rotated 0 times, done
	{
		if (rotlist_id == MAX_NUM_OF_ROTATIONS)
			return 1;

		//putting the num_of_threads_per_block pieces of most important rotations to the list
		for (parallel_rot_id=0; parallel_rot_id<num_of_threads_per_block; parallel_rot_id++)
		{
			if (number_of_req_rotations[parallel_rot_id] == 0)	//if the atom has not to be rotated anymore, dummy rotation
				new_rotlist_element = RLIST_DUMMY_MASK;
			else
			{
				atom_id = atom_id_of_numrots[parallel_rot_id];
				new_rotlist_element = ((int) atom_id) & RLIST_ATOMID_MASK;

				if (number_of_req_rotations[parallel_rot_id] == 1)
					new_rotlist_element |= RLIST_GENROT_MASK;
				else
				{
					rotbond_found = 0;
					rotbond_candidate = myligand->num_of_rotbonds - 1;
					remaining_rots_around_rotbonds = number_of_req_rotations[parallel_rot_id] - 1;	//-1 because of genrot

					while (rotbond_found == 0)
					{
						if (myligand->atom_rotbonds[atom_id][rotbond_candidate] != 0)	//if the atom has to be rotated around current candidate
						{
							if (remaining_rots_around_rotbonds == 1)	//if current value of remaining rots is 1, the proper rotbond is found
								rotbond_found = 1;
							else
								remaining_rots_around_rotbonds--;	//if not, decresing remaining rots (that is, skipping rotations which have to be performed later
						}

						if (rotbond_found == 0)
							rotbond_candidate--;

						if (rotbond_candidate < 0)
							return 1;
					}

					new_rotlist_element |= (((int) rotbond_candidate) << RLIST_RBONDID_SHIFT) & RLIST_RBONDID_MASK;
				}

				if (atom_wasnt_rotated_yet[atom_id] != 0)
					new_rotlist_element |= RLIST_FIRSTROT_MASK;

				//put atom_id's next rotation to rotlist
				atom_wasnt_rotated_yet[atom_id] = 0;
				(number_of_req_rotations[parallel_rot_id])--;
			}


			rotlist[rotlist_id] = new_rotlist_element;

			rotlist_id++;
		}

		make_reqrot_ordering(number_of_req_rotations, atom_id_of_numrots, myligand->num_of_atoms);
		(myligand->num_of_rotcyc)++;
	}

	// Subrotation lists
	// At some point, some code was proposed to build subrotations lists.
	// This was removed thereafter due to lower benefits in performance
	// compared to outer-loop pushing.
	// The printfs below are just vestige from that proposal,
	// which might be retaken eventually.
/*
	printf("\n# rotlist elements: %u\n", rotlist_id);
	for (unsigned int rot_cnt = 0; rot_cnt < myligand->num_of_rotations_required; rot_cnt++) {
		unsigned int atom_id = rotlist[rot_cnt] & RLIST_ATOMID_MASK;
		printf("rot-id: %u \tatom-id: %u\n", rot_cnt, atom_id);
	}

	printf("\n# atoms: %u\n", myligand->num_of_atoms);
	for (unsigned int atom_cnt = 0; atom_cnt < myligand->num_of_atoms; atom_cnt++) {
		printf("atom-id: %u \tnum-rot-req: %u\n", atom_cnt, number_of_req_rotations_copy[atom_cnt]);
	}
*/
	return 0;
}
