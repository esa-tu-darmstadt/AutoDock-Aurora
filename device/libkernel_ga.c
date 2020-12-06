#include "auxiliary.c"
#include "lga.c"

/*
IC:  initial calculation of energy of populations
GG:  genetic generation
LS:  local search
*/

// --------------------------------------------------------------------------
// Lamarckian Genetic-Algorithm (GA): GA + LS (Local Search)
// --------------------------------------------------------------------------
uint64_t libkernel_ga (
	const 	float*		PopulationCurrentInitial,
			uint64_t  	VEVMA_PopulationCurrentFinal,
			uint64_t  	VEVMA_EnergyCurrent,
			uint64_t	VEVMA_Evals_performed,
            uint64_t	VEVMA_Gens_performed,
			uint*    	dockpars_prng_states,
			uint		DockConst_pop_size,
		    uint        DockConst_num_of_energy_evals,
			uint        DockConst_num_of_generations,
		    float       DockConst_tournament_rate,
			float       DockConst_mutation_rate,
		    float       DockConst_abs_max_dmov,
			float       DockConst_abs_max_dang,
		    float       Host_two_absmaxdmov,
			float       Host_two_absmaxdang,
			float       DockConst_crossover_rate,
			uint        DockConst_num_of_lsentities,
			uchar       DockConst_num_of_genes,
	// pc
	const	int* 		PC_rotlist,
	const	int* 		PC_subrotlist_1,
	const	int* 		PC_subrotlist_2,
	const	int* 		PC_subrotlist_3,
	const	int* 		PC_subrotlist_4,
	const	int* 		PC_subrotlist_5,
	const	int* 		PC_subrotlist_6,
	const	int* 		PC_subrotlist_7,
	const	int* 		PC_subrotlist_8,
	const	int* 		PC_subrotlist_9,
	const	float*		PC_ref_coords_x,// TODO: merge them into a single one?
	const	float*		PC_ref_coords_y,
	const	float*		PC_ref_coords_z,
	const	float*		PC_rotbonds_moving_vectors,
	const	float*		PC_rotbonds_unit_vectors,
	const	float*		PC_ref_orientation_quats,
			uint		DockConst_rotbondlist_length,
			uint		subrotlist_1_length,
			uint		subrotlist_2_length,
			uint		subrotlist_3_length,
			uint		subrotlist_4_length,
			uint		subrotlist_5_length,
			uint		subrotlist_6_length,
			uint		subrotlist_7_length,
			uint		subrotlist_8_length,
			uint		subrotlist_9_length,
	// ia
	const 	float*		IA_IE_atom_charges,
	const	int*		IA_IE_atom_types,
	const	int*		IA_intraE_contributors,
	const	float*		IA_reqm,
	const	float*		IA_reqm_hbond,
	const	uint*		IA_atom1_types_reqm,
	const	uint*		IA_atom2_types_reqm,
	const	float*		IA_VWpars_AC,
	const	float*		IA_VWpars_BD,
	const	float*		IA_dspars_S,
	const	float*		IA_dspars_V,
			float		DockConst_smooth,
			uint		DockConst_num_of_intraE_contributors,
			float		DockConst_grid_spacing,
			uint		DockConst_num_of_atypes,
			float		DockConst_coeff_elec,
			float		DockConst_qasp,
			float		DockConst_coeff_desolv,
	// ie
	const	uint64_t	VEVMA_Fgrids,
			uchar		DockConst_g1,
			uint		DockConst_g2,
			uint		DockConst_g3,
			uchar		DockConst_num_of_atoms,
			float		DockConst_gridsize_x_minus1,
			float		DockConst_gridsize_y_minus1,
			float		DockConst_gridsize_z_minus1,
			uint		Host_mul_tmp2,
			uint		Host_mul_tmp3,
	// ls
			ushort		DockConst_max_num_of_iters,
			float		DockConst_rho_lower_bound,
			float		DockConst_base_dmov_mul_sqrt3,
			float		DockConst_base_dang_mul_sqrt3,
			uchar		DockConst_cons_limit,
	// Values changing every LGA run
            uint		Host_num_of_runs
)
{

	#pragma omp parallel for
	for (unsigned int run_cnt = 0; run_cnt < Host_num_of_runs; run_cnt++) {

/*
		printf(" %u", run_cnt+1); 
		fflush(stdout);
*/
		// Values changing every LGA run
		unsigned int uint_run_cnt  = run_cnt;
		unsigned int Host_Offset_Pop = run_cnt * DockConst_pop_size* ACTUAL_GENOTYPE_LENGTH;
		unsigned int Host_Offset_Ene = run_cnt * DockConst_pop_size;

        lga(
			PopulationCurrentInitial,
			VEVMA_PopulationCurrentFinal,
			VEVMA_EnergyCurrent,
			VEVMA_Evals_performed,
			VEVMA_Gens_performed,
			dockpars_prng_states,
			DockConst_pop_size,
			DockConst_num_of_energy_evals,
			DockConst_num_of_generations,
			DockConst_tournament_rate,
			DockConst_mutation_rate,
			DockConst_abs_max_dmov,
			DockConst_abs_max_dang,
			Host_two_absmaxdmov,
			Host_two_absmaxdang,
			DockConst_crossover_rate,
			DockConst_num_of_lsentities,
			DockConst_num_of_genes,
			// pc
			PC_rotlist,
			PC_subrotlist_1,
			PC_subrotlist_2,
			PC_subrotlist_3,
			PC_subrotlist_4,
			PC_subrotlist_5,
			PC_subrotlist_6,
			PC_subrotlist_7,
			PC_subrotlist_8,
			PC_subrotlist_9,
			PC_ref_coords_x,// TODO: merge them into a single one?
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			DockConst_rotbondlist_length,
			subrotlist_1_length,
			subrotlist_2_length,
			subrotlist_3_length,
			subrotlist_4_length,
			subrotlist_5_length,
			subrotlist_6_length,
			subrotlist_7_length,
			subrotlist_8_length,
			subrotlist_9_length,
			// ia
			IA_IE_atom_charges,
			IA_IE_atom_types,
			IA_intraE_contributors,
			IA_reqm,
			IA_reqm_hbond,
			IA_atom1_types_reqm,
			IA_atom2_types_reqm,
			IA_VWpars_AC,
			IA_VWpars_BD,
			IA_dspars_S,
			IA_dspars_V,
			DockConst_smooth,
			DockConst_num_of_intraE_contributors,
			DockConst_grid_spacing,
			DockConst_num_of_atypes,
			DockConst_coeff_elec,
			DockConst_qasp,
			DockConst_coeff_desolv,
			// ie
			VEVMA_Fgrids,
			DockConst_g1,
			DockConst_g2,
			DockConst_g3,
			DockConst_num_of_atoms,
			DockConst_gridsize_x_minus1,
			DockConst_gridsize_y_minus1,
			DockConst_gridsize_z_minus1,
			Host_mul_tmp2,
			Host_mul_tmp3,
			// ls
			DockConst_max_num_of_iters,
			DockConst_rho_lower_bound,
			DockConst_base_dmov_mul_sqrt3,
			DockConst_base_dang_mul_sqrt3,
			DockConst_cons_limit,
			// Values changing every LGA run
			uint_run_cnt,
			Host_Offset_Pop,
			Host_Offset_Ene
		);
		
	} // End of for loop

	return 0;
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
