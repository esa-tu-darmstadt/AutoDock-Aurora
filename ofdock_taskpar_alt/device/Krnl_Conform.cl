#include "auxiliary_conform.cl"

// --------------------------------------------------------------------------
// The function changes the conformation of myligand according to 
// the genotype given by the second parameter.
// Originally from: processligand.c
// --------------------------------------------------------------------------
__kernel __attribute__ ((reqd_work_group_size(1,1,1)))
void Krnl_Conform(__global const int*             restrict GlobLigand_rotbonds,
		  __global const int*             restrict GlobLigand_atom_rotbonds,
		  __global const float*           restrict GlobLigand_rotbonds_moving_vectors,
		  __global const float*           restrict GlobLigand_rotbonds_unit_vectors,
		  __global const Ligandconstant*  restrict LigConst)
{
#ifdef EMULATOR
	printf("Krnl Conform!!\n");
#endif
	// --------------------------------------------------------------
	// Wait for ligand data and genotypes
	// --------------------------------------------------------------
	__local float myligand_atom_idxyzq[MAX_NUM_OF_ATOMS*5];
	__local float genotype[40];

	uint init_cnt;

	for (init_cnt=0; init_cnt<LigConst->num_of_atoms*5; init_cnt++)
	{
		myligand_atom_idxyzq[init_cnt] = read_channel_altera(chan_GA2Conf_ligandatom_idxyzq);
	}

	for (init_cnt=0; init_cnt<40; init_cnt++)
	{
		genotype[init_cnt] = read_channel_altera(chan_GA2Conf_genotype);
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
	
	//moving ligand to origo
	get_movvec_to_origo(LigConst->num_of_atoms, myligand_atom_idxyzq, movvec_to_origo);
	move_ligand        (LigConst->num_of_atoms, myligand_atom_idxyzq, movvec_to_origo);

	//for each atom of the ligand
	// **********************************************
	// ADD VENDOR SPECIFIC PRAGMA
	// **********************************************
	LOOP_CHANGE_CONFORM_1:
	for (atom_id=0; atom_id<LigConst->num_of_atoms; atom_id++)						
	{
		#if defined (DEBUG_KERNEL_CHANGE_CONFORM)
		printf("\n\n\nROTATING atom %d ", atom_id);
		#endif

		//if the ligand has rotatable bonds
		if (LigConst->num_of_rotbonds != 0)
		{
			//for each rotatable bond
			// **********************************************
			// ADD VENDOR SPECIFIC PRAGMA
			// **********************************************	
			LOOP_CHANGE_CONFORM_2:
			for (rotbond_id=0; rotbond_id<LigConst->num_of_rotbonds; rotbond_id++)
			{
				//if the atom has to be rotated around this bond
				if (GlobLigand_atom_rotbonds [atom_id*32 + rotbond_id] != 0)
				{
					#if defined (DEBUG_KERNEL_CHANGE_CONFORM)
					printf("around rotatable bond %d\n", rotbond_id);
					#endif

					//rotating
					rotate_custom(&myligand_atom_idxyzq [atom_id*5+1], 
						      GlobLigand_rotbonds_moving_vectors [rotbond_id*3],
						      GlobLigand_rotbonds_moving_vectors [rotbond_id*3+1],
						      GlobLigand_rotbonds_moving_vectors [rotbond_id*3+2],
						      GlobLigand_rotbonds_unit_vectors   [rotbond_id*3],
						      GlobLigand_rotbonds_unit_vectors   [rotbond_id*3+1],
						      GlobLigand_rotbonds_unit_vectors   [rotbond_id*3+2],
						      genotype [6+rotbond_id]);

				} // End of if (GlobLigand_atom_rotbonds [atom_id*32 + rotbond_id] != 0)
			} // End of for-loop (rotbond_id)	
		} // End of if (myligand_num_of_rotbonds != 0)

		#if defined (DEBUG_KERNEL_CHANGE_CONFORM)
		printf("according to general rotation\n");
		#endif

		//general rotation
		rotate_custom(&myligand_atom_idxyzq [atom_id*5+1], 
		              genrot_movvec[0],
			      genrot_movvec[1],
			      genrot_movvec[2],
		              genrot_unitvec[0],
			      genrot_unitvec[1],
			      genrot_unitvec[2],
		              genotype [5]);

	} // End of for-loop (atom_id)

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

	// --------------------------------------------------------------
	// Send ligand data to InterE and IntraE Kernels
	// --------------------------------------------------------------
	for (init_cnt=0; init_cnt<LigConst->num_of_atoms*5; init_cnt++)
	{
		write_channel_altera(chan_Conf2Intere_ligandatom_idxyzq, myligand_atom_idxyzq[init_cnt]);
		write_channel_altera(chan_Conf2Intrae_ligandatom_idxyzq, myligand_atom_idxyzq[init_cnt]);
	}
	// --------------------------------------------------------------

}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
