// --------------------------------------------------------------------------
// Originally from: processligand.c
// --------------------------------------------------------------------------
__kernel
void Krnl_IntraE(
             __global const float*           restrict GlobFgrids,
	     __global       float*           restrict GlobPopulationCurrent,
	     __global       float*           restrict GlobEnergyCurrent,
	     __global       float*           restrict GlobPopulationNext,
	     __global       float*           restrict GlobEnergyNext,
             __global       unsigned int*    restrict GlobPRNG,
	     __global const kernelconstant*  restrict KerConst,
	     __global const Dockparameters*  restrict DockConst)
{

	__local float loc_coords_x[MAX_NUM_OF_ATOMS];
	__local float loc_coords_y[MAX_NUM_OF_ATOMS];
	__local float loc_coords_z[MAX_NUM_OF_ATOMS];

	char active = 1;
	char mode   = 0;
	uint cnt    = 0;   

	int contributor_counter;
	char atom1_id, atom2_id, atom1_typeid, atom2_typeid;
	float subx, suby, subz, distance_leo;

 	// Altera doesn't support power function 	
	// so this is implemented with multiplications 	
	// Full precision is used 	
	float distance_pow_2, distance_pow_4, distance_pow_6, distance_pow_10, distance_pow_12;
	float intraE;

while(active) {
	//printf("BEFORE In INTRA CHANNEL\n");
	// --------------------------------------------------------------
	// Wait for ligand atomic coordinates in channel
	// --------------------------------------------------------------

	for (uint pipe_cnt=0; pipe_cnt<DockConst->num_of_atoms; pipe_cnt++) {
		loc_coords_x[pipe_cnt] = read_channel_altera(chan_Conf2Intrae_x);
		mem_fence(CLK_CHANNEL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		loc_coords_y[pipe_cnt] = read_channel_altera(chan_Conf2Intrae_y);
		mem_fence(CLK_CHANNEL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		loc_coords_z[pipe_cnt] = read_channel_altera(chan_Conf2Intrae_z);
		mem_fence(CLK_CHANNEL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		active = read_channel_altera(chan_Conf2Intrae_active);
		mem_fence(CLK_CHANNEL_MEM_FENCE);
		mode   = read_channel_altera(chan_Conf2Intrae_mode);
		mem_fence(CLK_CHANNEL_MEM_FENCE);
		cnt    = read_channel_altera(chan_Conf2Intrae_cnt);
	}
	// --------------------------------------------------------------
	//printf("AFTER In INTRA CHANNEL\n");

	if (active == 0) {printf("	%-20s: %s\n", "Krnl_IntraE", "disabled");}

	intraE = 0.0f;

	//for each intramolecular atom contributor pair
	// **********************************************
	// ADD VENDOR SPECIFIC PRAGMA
	// **********************************************
	LOOP_INTRAE_1:
	for (uint contributor_counter=0; contributor_counter<DockConst->num_of_intraE_contributors; contributor_counter++)
	{	
		atom1_id = KerConst->intraE_contributors_const[3*contributor_counter]; 
		atom2_id = KerConst->intraE_contributors_const[3*contributor_counter+1];

		subx = loc_coords_x[atom1_id] - loc_coords_x[atom2_id];
		suby = loc_coords_y[atom1_id] - loc_coords_y[atom2_id];
		subz = loc_coords_z[atom1_id] - loc_coords_z[atom2_id];
		distance_leo = sqrt(subx*subx + suby*suby + subz*subz)*DockConst->grid_spacing;

		if (distance_leo < 1.0f) {
			#if defined (DEBUG_KERNEL_INTRA_E)
			printf("\n\nToo low distance (%f) between atoms %u and %u\n", distance_leo, atom1_id, atom2_id);
			#endif
			//return HIGHEST_ENERGY;	//returning maximal value
			distance_leo = 1.0f;
		}

		#if defined (DEBUG_KERNEL_INTRA_E)
		printf("\n\nCalculating energy contribution of atoms %u and %u\n", atom1_id+1, atom2_id+1);
		printf("Distance: %f\n", distance_leo);
		#endif

		distance_pow_2  = distance_leo*distance_leo; 		
		distance_pow_4  = distance_pow_2*distance_pow_2; 		
		distance_pow_6  = distance_pow_2*distance_pow_4; 		
		distance_pow_10 = distance_pow_4*distance_pow_6; 		
		distance_pow_12 = distance_pow_6*distance_pow_6;
		
		//but only if the distance is less than distance cutoff value and 20.48A (because of the tables)
		if ((distance_leo < 8.0f) && (distance_leo < 20.48f)) 
		{
			atom1_typeid = KerConst->atom_types_const [atom1_id];
			atom2_typeid = KerConst->atom_types_const [atom2_id];

			//calculating van der Waals / hydrogen bond term
			intraE += KerConst->VWpars_AC_const[atom1_typeid * DockConst->num_of_atypes+atom2_typeid]/distance_pow_12;

			if (KerConst->intraE_contributors_const[3*contributor_counter+2] == 1)	//H-bond
				intraE-= KerConst->VWpars_BD_const[atom1_typeid*DockConst->num_of_atypes+atom2_typeid]/distance_pow_10;	
			else	//van der Waals
				intraE-= KerConst->VWpars_BD_const[atom1_typeid*DockConst->num_of_atypes+atom2_typeid]/distance_pow_6;

			//calculating electrostatic term
			intraE+= DockConst->coeff_elec*KerConst->atom_charges_const[atom1_id]*KerConst->atom_charges_const[atom2_id]/(distance_leo*(-8.5525f + 86.9525f/(1.0f + 7.7839f*exp(-0.3154f*distance_leo))));

			//calculating desolvation term
			intraE+= (
				  ( KerConst->dspars_S_const[atom1_typeid] + DockConst->qasp*fabs(KerConst->atom_charges_const[atom1_id]) ) * KerConst->dspars_V_const[atom2_typeid] + 
				  ( KerConst->dspars_S_const[atom2_typeid] + DockConst->qasp*fabs(KerConst->atom_charges_const[atom2_id]) ) * KerConst->dspars_V_const[atom1_typeid]) * 
				 DockConst->coeff_desolv*exp(-distance_leo*distance_leo/25.92f);
	
		} // End of if: if ((dist < dcutoff) && (dist < 20.48))	

	} // End of LOOP_INTRAE_1

	// --------------------------------------------------------------
	// Send intramolecular energy to channel
	// --------------------------------------------------------------
	write_channel_altera(chan_Intrae2Store_intrae, intraE);
	mem_fence(CLK_CHANNEL_MEM_FENCE);
	write_channel_altera(chan_Intrae2Store_active, active);
	mem_fence(CLK_CHANNEL_MEM_FENCE);
	write_channel_altera(chan_Intrae2Store_mode,   mode);
	mem_fence(CLK_CHANNEL_MEM_FENCE);
	write_channel_altera(chan_Intrae2Store_cnt,    cnt);
	// --------------------------------------------------------------

	} // End of while(1)
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
