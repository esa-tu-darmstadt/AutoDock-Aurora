// --------------------------------------------------------------------------
// Perform Genetic Generation (GG)
// --------------------------------------------------------------------------

__kernel __attribute__ ((max_global_work_dim(0)))
void Krnl_GG(__global const float*           restrict GlobPopulationCurrent,
	     __global const float*           restrict GlobEnergyCurrent,
	     __global       float*           restrict GlobPopulationNext,
	     __global       float*           restrict GlobEnergyNext,
             __global       unsigned int*    restrict GlobPRNG,
	     __constant     Dockparameters*  restrict DockConst      			      	      
)
{		
	      char active = 1; 	
	const char mode   = 2; 	
	      char ack    = 0;
	
	// Find_best 	
	uint best_entity_id; 	
	__local float loc_energies[MAX_POPSIZE]; 	
	
	// Binary tournament 	
	uint parent1, parent2; 

	/*__local*/ float local_entity_1     [ACTUAL_GENOTYPE_LENGTH]; 	
	/*__local*/ float local_entity_2     [ACTUAL_GENOTYPE_LENGTH];	
/*
	float __attribute__((register)) local_entity_1 [ACTUAL_GENOTYPE_LENGTH]; 	
	float __attribute__((register)) local_entity_2 [ACTUAL_GENOTYPE_LENGTH];	
*/	
	
	__local float offspring_genotype [ACTUAL_GENOTYPE_LENGTH]; 	
	
while(active) {
	active = read_channel_altera(chan_GA2GG_active);

	#if defined (DEBUG_ACTIVE_KERNEL)
	if (active == 0) {printf("	%-20s: %s\n", "Krnl_GG", "must be disabled");}
	#endif

	for (ushort i=0; i<DockConst->pop_size; i++) {
		loc_energies[i] = GlobEnergyCurrent[i];
	}

	//Identifying best entity 		
	best_entity_id = find_best(loc_energies, DockConst->pop_size);

	//elitism - copying the best entity to new population 		
	for (uchar i=0; i<ACTUAL_GENOTYPE_LENGTH; i++) { 			
		GlobPopulationNext[i] = GlobPopulationCurrent[best_entity_id*ACTUAL_GENOTYPE_LENGTH+i]; 		
	} 		
	GlobEnergyNext[0] = loc_energies[best_entity_id];

	//new population consists of one member currently 		
	//new_pop_cnt = 1;

	//read GlobPRNG
	uint prng = GlobPRNG[0];

	for (ushort new_pop_cnt = 1; new_pop_cnt < DockConst->pop_size; new_pop_cnt++) {

		//selecting two individuals randomly 			
		binary_tournament_selection(&prng, loc_energies, &parent1, &parent2,			    
				            DockConst->pop_size, DockConst->tournament_rate);

		//mating parents				
		for (uchar i=0; i<ACTUAL_GENOTYPE_LENGTH; i++) {
			local_entity_1[i] = GlobPopulationCurrent[parent1*ACTUAL_GENOTYPE_LENGTH+i];
			local_entity_2[i] = GlobPopulationCurrent[parent2*ACTUAL_GENOTYPE_LENGTH+i];
		}
		
		// first two args are population [parent1], population [parent2] 			
		gen_new_genotype(&prng, local_entity_1, local_entity_2, 					 
				 DockConst->mutation_rate, DockConst->abs_max_dmov, DockConst->abs_max_dang,			 
				 DockConst->crossover_rate, offspring_genotype);		

		//////======================================================
		/*
		if ( (active!=0) || ((active==0) && (new_pop_cnt==1)) ) {

			write_channel_altera(chan_GG2Conf_active, active);
			mem_fence(CLK_CHANNEL_MEM_FENCE);
			write_channel_altera(chan_GG2Conf_mode,   mode);
			mem_fence(CLK_CHANNEL_MEM_FENCE);
			write_channel_altera(chan_GG2Conf_cnt,    new_pop_cnt);
			mem_fence(CLK_CHANNEL_MEM_FENCE);

			for (uchar pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
			GlobPopulationNext [new_pop_cnt*ACTUAL_GENOTYPE_LENGTH + pipe_cnt] = offspring_genotype [pipe_cnt];
			write_channel_altera(chan_GG2Conf_genotype, offspring_genotype[pipe_cnt]);
			}

		} // End of if enclosing write to channels
		*/
		//////======================================================
			
		write_channel_altera(chan_GG2Conf_active, active);
		mem_fence(CLK_CHANNEL_MEM_FENCE);
		write_channel_altera(chan_GG2Conf_mode,   mode);
		mem_fence(CLK_CHANNEL_MEM_FENCE);
		write_channel_altera(chan_GG2Conf_cnt,    new_pop_cnt);
		mem_fence(CLK_CHANNEL_MEM_FENCE);

		for (uchar pipe_cnt=0; pipe_cnt<ACTUAL_GENOTYPE_LENGTH; pipe_cnt++) {
		GlobPopulationNext [new_pop_cnt*ACTUAL_GENOTYPE_LENGTH + pipe_cnt] = offspring_genotype [pipe_cnt];
		write_channel_altera(chan_GG2Conf_genotype, offspring_genotype[pipe_cnt]);
		}

	} // End of for-loop new_pop_cnt

	//write back to GlobPRNG
	GlobPRNG[0] = prng;

	ack = read_channel_altera(chan_Store2GG_ack);

	mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_CHANNEL_MEM_FENCE); // lvs added during hw evaluation

	write_channel_altera(chan_GG2GA_eval_cnt, (DockConst->pop_size - 1));

} // End of while (active)

#if defined (DEBUG_ACTIVE_KERNEL)
printf("	%-20s: %s\n", "Krnl_GG", "disabled");
#endif

}
