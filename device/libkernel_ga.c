#include "auxiliary.h"
#include "rand_gen.h"
#include "lga.h"

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
	const	float*		PC_ref_coords_x,// TODO: merge them into a single one?
	const	float*		PC_ref_coords_y,
	const	float*		PC_ref_coords_z,
	const	float*		PC_rotbonds_moving_vectors,
	const	float*		PC_rotbonds_unit_vectors,
	const	float*		PC_ref_orientation_quats,
        uint		DockConst_rotbondlist_length,
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
	randf_vec_init();

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
			PC_ref_coords_x,// TODO: merge them into a single one?
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			DockConst_rotbondlist_length,
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

	randf_vec_fini();
	return 0;
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
