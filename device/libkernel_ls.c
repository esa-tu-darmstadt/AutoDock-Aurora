#include "defines.h"
#include "math.h"
#include "auxiliary.c"

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
void Krnl_LS(
		unsigned short            DockConst_max_num_of_iters,
		float                     DockConst_rho_lower_bound,
		float                     DockConst_base_dmov_mul_sqrt3,
		unsigned char             DockConst_num_of_genes,
   		float                     DockConst_base_dang_mul_sqrt3,
		unsigned char             DockConst_cons_limit
)
{	
	#if 0
	printf("\nLS: DockConst_max_num_of_iters: %u\n",   DockConst_max_num_of_iters);
	printf("LS: DockConst_rho_lower_bound: %f\n",      DockConst_rho_lower_bound);
	printf("LS: DockConst_base_dmov_mul_sqrt3: %f\n",  DockConst_base_dmov_mul_sqrt3);
	printf("LS: DockConst_num_of_genes: %u\n",  	   DockConst_num_of_genes);
	printf("LS: DockConst_base_dang_mul_sqrt3: %f\n",  DockConst_base_dang_mul_sqrt3);
	printf("LS: DockConst_cons_limit: %u\n",           DockConst_cons_limit);
	#endif

	float current_energy;
	float   genotype [ACTUAL_GENOTYPE_LENGTH];

		// LOOP_FOR_LS_READ_INPUT_GENOTYPE
		for (unsigned char i=0; i<DockConst_num_of_genes; i++) {
			read_pipe_block(pipe00ga2ls00ls100genotype, &genotype [i]);
		}
	
		float rho = 1.0f;
		unsigned short iteration_cnt = 0;
		unsigned char  cons_succ     = 0;
		unsigned char  cons_fail     = 0;
		unsigned int   LS_eval       = 0;
		boolean positive_direction = True;

		// performing local search
		// LOOP_WHILE_LS_ITERATION_RHO
		while ((iteration_cnt < DockConst_max_num_of_iters) && (rho > DockConst_rho_lower_bound)) {
			// -----------------------------------------------
			// Exit condition is groups here. It allows pipelining
			if (positive_direction == True) { 
				if (cons_succ >= DockConst_cons_limit) {
					rho = LS_EXP_FACTOR*rho;
					cons_fail = 0;
					cons_succ = 0;
				}
				else if (cons_fail >= DockConst_cons_limit) {
					rho = LS_CONT_FACTOR*rho;
					cons_fail = 0;
					cons_succ = 0;
				}
				iteration_cnt++;
			}

			#if defined (DEBUG_KRNL_LS1)
			printf("LS1 positive?: %u, iteration_cnt: %u, rho: %f, limit rho: %f\n", positive_direction, iteration_cnt, rho, DockConst_rho_lower_bound);
			#endif
			// -----------------------------------------------

			float entity_possible_new_genotype   [ACTUAL_GENOTYPE_LENGTH];
			float genotype_bias                  [ACTUAL_GENOTYPE_LENGTH];
			float deviate_plus_bias              [ACTUAL_GENOTYPE_LENGTH];
			float deviate_minus_bias             [ACTUAL_GENOTYPE_LENGTH];

			// Tell Krnl_Conf_Arbiter, LS1 is done
			// Not completely strict as the (iteration_cnt < DockConst_max_num_of_iters) is ignored
			// In practice, rho condition dominates most of the cases
			int tmp_int = (rho < DockConst_rho_lower_bound)?0:1;
			write_pipe_block(pipe00ls2arbiter00ls100end, &tmp_int);

			// new random deviate
			// rho is the deviation of the uniform distribution
			// LOOP_FOR_LS_WRITE_GENOTYPE
			for (unsigned char i=0; i<DockConst_num_of_genes; i++) {
				float tmp_prng;
				read_pipe_block(pipe00prng2ls00float00prng, &tmp_prng);

				// tmp1 is genotype_deviate
				float tmp1 = rho * (2.0f*tmp_prng - 1.0f);

				if (i<3) { tmp1 = tmp1 * DockConst_base_dmov_mul_sqrt3; }
				else     { tmp1 = tmp1 * DockConst_base_dang_mul_sqrt3; }

				float deviate = 0.4f*tmp1;

				// tmp2 is the addition: genotype_deviate + genotype_bias
				// tmp3 is entity_possible_new_genotype
				float tmp_bias = (iteration_cnt == 1)? 0.0f:genotype_bias[i];
				float bias = 0.6f * tmp_bias;

				deviate_plus_bias  [i] = deviate + bias;
				deviate_minus_bias [i] = deviate - bias;

				float tmp2 = tmp1 + tmp_bias;
				float tmp3 = (positive_direction == True)? (genotype [i] + tmp2): (genotype [i] - tmp2);

				if (i>2) {if (i==4) { tmp3 = map_angle_180(tmp3); }
					  else      { tmp3 = map_angle_360(tmp3); }}

				entity_possible_new_genotype [i] = tmp3;
				write_pipe_block(pipe00ls2conf00ls100genotype, &tmp3);

				#if defined (DEBUG_KRNL_LS1)
				printf("LS1_genotype sent\n");
				#endif
			}

			//printf("Energy to calculate sent from LS ... ");

			float energyIA_LS_rx;
			float energyIE_LS_rx;

/*
			intra_valid = read_pipe(pipe00intrae2storels00ls100intrae, &energyIA_LS_rx);
			inter_valid = read_pipe(pipe00intere2storels00ls100intere, &energyIE_LS_rx);
*/

			float candidate_energy = energyIA_LS_rx + energyIE_LS_rx;

			// update LS energy-evaluation count
			LS_eval++;

			if (candidate_energy < current_energy) {
				// updating offspring_genotype
				// updating genotype_bias

				// LOOP_FOR_LS_FLOATPT_UPDATE_POS_GENOTYPE
				for (unsigned char i=0; i<DockConst_num_of_genes; i++) {
					genotype_bias [i] = (positive_direction == True) ? deviate_plus_bias  [i] : 
											   deviate_minus_bias [i] ;
					genotype [i] = entity_possible_new_genotype [i];
				}	

				current_energy = candidate_energy;
				cons_succ++;
				cons_fail = 0;
				positive_direction = True;				
			}
			else {
				// updating (halving) genotype_bias

				// LOOP_FOR_LS_FLOATPT_UPDATE_NEG_GENOTYPE
				for (unsigned char i=0; i<DockConst_num_of_genes; i++) {
					genotype_bias [i] = (iteration_cnt == 1)? 0.0f: (0.5f*genotype_bias [i]);
				}

				if (positive_direction == False) {
					cons_fail++;
					cons_succ = 0;
				}
				positive_direction = (positive_direction == True) ? False: True;
			}

		} // end of while (iteration_cnt) && (rho)

		#if defined (DEBUG_KRNL_LS1)
		printf("Out of while iter LS1\n");
		#endif
		
		// write back data to GA
		// LOOP_FOR_LS_WRITEBACK2GA
		for (unsigned char i=0; i<DockConst_num_of_genes; i++) {
			if (i == 0) {
				float2 evalenergy  = {*(float*)&LS_eval, current_energy};
				write_pipe_block(pipe00ls2ga00ls100evalenergy, &evalenergy);	
			}
			write_pipe_block(pipe00ls2ga00ls100genotype, &genotype [i]);
		}

#if defined (DEBUG_ACTIVE_KERNEL)
printf("	%-20s: %s\n", "Krnl_LS1", "disabled");		
#endif
	
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
