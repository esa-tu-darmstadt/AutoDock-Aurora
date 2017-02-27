#include "perform_elitist_selection.cl"

__kernel 
void __attribute__ ((reqd_work_group_size(NUM_OF_THREADS_PER_BLOCK,1,1)))
gen_and_eval_newpops(      char   dockpars_num_of_atoms,
		           char   dockpars_num_of_atypes,
		           int    dockpars_num_of_intraE_contributors,
		           char   dockpars_gridsize_x,
       		           char   dockpars_gridsize_y,
        	           char   dockpars_gridsize_z,
        	           float  dockpars_grid_spacing,	
            __global const float* restrict dockpars_fgrids, // cannot be allocated in __constant (too large)
                           int    dockpars_rotbondlist_length,
  	                   float  dockpars_coeff_elec,
		           float  dockpars_coeff_desolv,
            __global const float* restrict dockpars_conformations_current,
                  __global float* restrict dockpars_energies_current,
                  __global float* restrict dockpars_conformations_next,
                  __global float* restrict dockpars_energies_next,
                  __global int*   restrict dockpars_evals_of_new_entities,
                  __global unsigned int* restrict dockpars_prng_states,
		           int    dockpars_pop_size,
		           int    dockpars_num_of_genes,
		           float  dockpars_tournament_rate,
	                   float  dockpars_crossover_rate,
		           float  dockpars_mutation_rate,
		           float  dockpars_abs_max_dmov,
		           float  dockpars_abs_max_dang,
		           float  dockpars_qasp,

__global kernelconstant* restrict KerConst

/*
	        __constant float* restrict atom_charges_const,
                __constant char*  restrict atom_types_const,
	        __constant char*  restrict intraE_contributors_const,
                __constant float* restrict VWpars_AC_const,
                __constant float* restrict VWpars_BD_const,
                __constant float* restrict dspars_S_const,
                __constant float* restrict dspars_V_const,
                __constant int*   restrict rotlist_const,
                __constant float* restrict ref_coords_x_const,
                __constant float* restrict ref_coords_y_const,
                __constant float* restrict ref_coords_z_const,
                __constant float* restrict rotbonds_moving_vectors_const,
                __constant float* restrict rotbonds_unit_vectors_const,
                __constant float* restrict ref_orientation_quats_const

*/

)
//The GPU global function
{

	__local float offspring_genotype[ACTUAL_GENOTYPE_LENGTH];
	__local int parent_candidates[4];
	__local float candidate_energies[4];
	__local int parents[2];
	__local int run_id;
	__local int covr_point[2];
	__local float randnums[10];
	int temp_covr_point;
	int gene_counter;
	__local float energy;	//could be shared since only thread 0 will use it

	// Altera doesn't allow local var outside kernels
	// so this local vars are passed from a kernel
	__local float best_energies[NUM_OF_THREADS_PER_BLOCK];
	__local int best_IDs[NUM_OF_THREADS_PER_BLOCK];
        __local int best_ID[1]; //__local int best_ID;

	// Altera doesn't allow local var outside kernels
	// so this local vars are passed from a kernel
	__local float calc_coords_x[MAX_NUM_OF_ATOMS];
	__local float calc_coords_y[MAX_NUM_OF_ATOMS];
	__local float calc_coords_z[MAX_NUM_OF_ATOMS];
	__local float partial_energies[NUM_OF_THREADS_PER_BLOCK];

	//in this case this block is responsible for elitist selection
	if ((get_group_id(0) % dockpars_pop_size) == 0)	

		//gpu_perform_elitist_selection(&dockpars);
		gpu_perform_elitist_selection(dockpars_pop_size,
					      dockpars_energies_current,
					      dockpars_energies_next,
					      dockpars_evals_of_new_entities,
					      dockpars_num_of_genes,
					      dockpars_conformations_next,
				              dockpars_conformations_current,
					      best_energies,
					      best_IDs,
					      best_ID);
	else
	{
		//generating the following random numbers: [0..3] for parent candidates, 
		//[4..5] for binary tournaments, [6] for deciding crossover,
		//[7..8] for crossover points, [9] for local search

		for (gene_counter=get_local_id(0);
		     gene_counter<10;
		     gene_counter+=NUM_OF_THREADS_PER_BLOCK)
			randnums[gene_counter] = gpu_randf(dockpars_prng_states);

		//determining run ID
		if (get_local_id(0) == 0)
			run_id = get_group_id(0) / dockpars_pop_size;

		//performing binary tournament selection
		barrier(CLK_LOCAL_MEM_FENCE);

		if (get_local_id(0) < 4)	//it is not ensured that the four candidates will be different...
		{
			//using randnums[0..3]
			parent_candidates[get_local_id(0)] = (int) (dockpars_pop_size*randnums[get_local_id(0)]);
			candidate_energies[get_local_id(0)] = dockpars_energies_current[run_id*dockpars_pop_size+
											parent_candidates[get_local_id(0)]];
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		if (get_local_id(0) < 2)
		{
			if (candidate_energies[2*get_local_id(0)] < candidate_energies[2*get_local_id(0)+1])
				//using randnum[4..5]
				if (100.0f*randnums[4+get_local_id(0)] < dockpars_tournament_rate)
					parents[get_local_id(0)] = parent_candidates[2*get_local_id(0)];
				else
					parents[get_local_id(0)] = parent_candidates[2*get_local_id(0)+1];
			else
				if (100.0f*randnums[4+get_local_id(0)] < dockpars_tournament_rate)
					parents[get_local_id(0)] = parent_candidates[2*get_local_id(0)+1];
				else
					parents[get_local_id(0)] = parent_candidates[2*get_local_id(0)];
		}

		//performing crossover
		barrier(CLK_LOCAL_MEM_FENCE);

		if (100.0f*randnums[6] < dockpars_crossover_rate)	//using randnums[6]
		{
			if (get_local_id(0) < 2)
				//using randnum[7..8]
				covr_point[get_local_id(0)] = (int) ((dockpars_num_of_genes-1)*randnums[7+get_local_id(0)]);

			barrier(CLK_LOCAL_MEM_FENCE);
			if (get_local_id(0) == 0)	//covr_point[0] should store the lower crossover-point
				if (covr_point[1] < covr_point[0])
				{
					temp_covr_point = covr_point[1];
					covr_point[1] = covr_point[0];
					covr_point[0] = temp_covr_point;
				}

			barrier(CLK_LOCAL_MEM_FENCE);

			for (gene_counter=get_local_id(0);
			     gene_counter<dockpars_num_of_genes;
			     gene_counter+=NUM_OF_THREADS_PER_BLOCK)
			{
				if (covr_point[0] != covr_point[1])	//two-point crossover
					if ((gene_counter <= covr_point[0]) || (gene_counter > covr_point[1]))
						offspring_genotype[gene_counter] = dockpars_conformations_current[(run_id*dockpars_pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter];
					else
						offspring_genotype[gene_counter] = dockpars_conformations_current[(run_id*dockpars_pop_size+parents[1])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter];
				else									//single-point crossover
					if (gene_counter <= covr_point[0])
						offspring_genotype[gene_counter] = dockpars_conformations_current[(run_id*dockpars_pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter];
					else
						offspring_genotype[gene_counter] = dockpars_conformations_current[(run_id*dockpars_pop_size+parents[1])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter];
			}

		}
		else	//no crossover
		{
#if defined (ASYNC_COPY)
			async_work_group_copy(offspring_genotype,
        	        	              dockpars_conformations_current+(run_id*dockpars_pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM,
                	        	      dockpars_num_of_genes,0);

#else
			for (gene_counter=get_local_id(0);
			     gene_counter<dockpars_num_of_genes;
			     gene_counter+=NUM_OF_THREADS_PER_BLOCK)
				offspring_genotype[gene_counter] = dockpars_conformations_current[(run_id*dockpars_pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter];
#endif
		}


		// L30nardoSV added
		// This was apparently missing
		barrier(CLK_LOCAL_MEM_FENCE);

		//performing mutation
		for (gene_counter=get_local_id(0);
		     gene_counter<dockpars_num_of_genes;
		     gene_counter+=NUM_OF_THREADS_PER_BLOCK)
		{		
			if (100.0f*gpu_randf(dockpars_prng_states) < dockpars_mutation_rate)
			{
				if (gene_counter < 3)
					offspring_genotype[gene_counter] += dockpars_abs_max_dmov*(2*gpu_randf(dockpars_prng_states)-1);
				else
				{
					offspring_genotype[gene_counter] += dockpars_abs_max_dang*(2*gpu_randf(dockpars_prng_states)-1);
					map_angle(&(offspring_genotype[gene_counter]));
				}
			}
		}

		//calculating energy of new offspring
		barrier(CLK_LOCAL_MEM_FENCE);


		// =================================================================
		// =================================================================
		// =================================================================
		//WARNING: only energy of threadIdx.x=0 will be valid
		//gpu_calc_energy(&dockpars, offspring_genotype, &energy, &run_id);
		calc_energy(    dockpars_rotbondlist_length,
				dockpars_num_of_atoms,
				dockpars_gridsize_x,
				dockpars_gridsize_y,
	                    	dockpars_gridsize_z,
				dockpars_fgrids,
				dockpars_num_of_atypes,
				dockpars_num_of_intraE_contributors,
				dockpars_grid_spacing,
				dockpars_coeff_elec,
                	    	dockpars_qasp,
				dockpars_coeff_desolv,
				offspring_genotype,
				&energy,
				&run_id,
   				// Altera doesn't allow local var outside kernels
				// so this local vars are passed from a kernel
	               		calc_coords_x,
	                	calc_coords_y,
	                	calc_coords_z,
	                	partial_energies,
				
				KerConst

/*
                        	atom_charges_const,
	                   	atom_types_const,
        	            	intraE_contributors_const,
                	    	VWpars_AC_const,
                        	VWpars_BD_const,
	                    	dspars_S_const,
        	            	dspars_V_const,
                	    	rotlist_const,
                        	ref_coords_x_const,
	                    	ref_coords_y_const,
        	            	ref_coords_z_const,
                	    	rotbonds_moving_vectors_const,
                        	rotbonds_unit_vectors_const,
	                    	ref_orientation_quats_const
*/
);
		// =================================================================
		// =================================================================
		// =================================================================


		if (get_local_id(0) == 0)
		{
			dockpars_evals_of_new_entities[get_group_id(0)] = 1;
		}


		if (get_local_id(0) == 0)
		{
			dockpars_energies_next[get_group_id(0)] = energy;
		}

		//copying new offspring to next generation

#if defined (ASYNC_COPY)
	        async_work_group_copy(dockpars_conformations_next + GENOTYPE_LENGTH_IN_GLOBMEM*get_group_id(0),
        	                      offspring_genotype,
                	              dockpars_num_of_genes,0);
#else
		for (gene_counter=get_local_id(0);
		     gene_counter<dockpars_num_of_genes;
		     gene_counter+=NUM_OF_THREADS_PER_BLOCK)
			dockpars_conformations_next[GENOTYPE_LENGTH_IN_GLOBMEM*get_group_id(0)+gene_counter] = offspring_genotype[gene_counter];
#endif
	}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
        if (get_local_id(0) == 0) {

		printf("%u,lid: %zu, groupid: %zu, #groups: %zu, dockpars_energies_next[%zu]: %f\n",
			GENOTYPE_LENGTH_IN_GLOBMEM,			
			get_local_id(0),get_group_id(0),get_num_groups(0),
			get_group_id(0),dockpars_energies_current[get_group_id(0)]);	
      	}
*/
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}
