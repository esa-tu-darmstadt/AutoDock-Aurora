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
	const 	uint64_t	VEVMA_PopulationCurrentInitial,
			uint64_t  	VEVMA_PopulationCurrentFinal,
			uint64_t  	VEVMA_EnergyCurrent,
			uint64_t	VEVMA_Evals_performed,
            uint64_t	VEVMA_Gens_performed,
			uint64_t    VEVMA_dockpars_prng_states,
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
	const	uint64_t	VEVMA_pc_rotlist,
	const	uint64_t	VEVMA_pc_ref_coords_x,// TODO: merge them into a single one?
	const	uint64_t	VEVMA_pc_ref_coords_y,
	const	uint64_t	VEVMA_pc_ref_coords_z,
	const	uint64_t	VEVMA_pc_rotbonds_moving_vectors,
	const	uint64_t	VEVMA_pc_rotbonds_unit_vectors,
	const	uint64_t	VEVMA_pc_ref_orientation_quats,
			uint		DockConst_rotbondlist_length,
	// ia
	const 	uint64_t	VEVMA_ia_ie_atom_charges,
	const	uint64_t	VEVMA_ia_ie_atom_types,
	const	uint64_t	VEVMA_ia_intraE_contributors,
	const	uint64_t	VEVMA_ia_reqm,
	const	uint64_t	VEVMA_ia_reqm_hbond,
	const	uint64_t	VEVMA_ia_atom1_types_reqm,
	const	uint64_t	VEVMA_ia_atom2_types_reqm,
	const	uint64_t	VEVMA_ia_VWpars_AC,
	const	uint64_t	VEVMA_ia_VWpars_BD,
	const	uint64_t	VEVMA_ia_dspars_S,
	const	uint64_t	VEVMA_ia_dspars_V,
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

	for (unsigned int run_cnt = 0; run_cnt < Host_num_of_runs; run_cnt++) {
		printf(" %u", run_cnt+1); 
		fflush(stdout);

		// Values changing every LGA run
		unsigned int uint_run_cnt  = run_cnt;
		unsigned int Host_Offset_Pop = run_cnt * DockConst_pop_size* ACTUAL_GENOTYPE_LENGTH;
		unsigned int Host_Offset_Ene = run_cnt * DockConst_pop_size;

        lga(
			VEVMA_PopulationCurrentInitial,
			VEVMA_PopulationCurrentFinal,
			VEVMA_EnergyCurrent,
			VEVMA_Evals_performed,
			VEVMA_Gens_performed,
			VEVMA_dockpars_prng_states,
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
			VEVMA_pc_rotlist,
			VEVMA_pc_ref_coords_x,// TODO: merge them into a single one?
			VEVMA_pc_ref_coords_y,
			VEVMA_pc_ref_coords_z,
			VEVMA_pc_rotbonds_moving_vectors,
			VEVMA_pc_rotbonds_unit_vectors,
			VEVMA_pc_ref_orientation_quats,
			DockConst_rotbondlist_length,
			// ia
			VEVMA_ia_ie_atom_charges,
			VEVMA_ia_ie_atom_types,
			VEVMA_ia_intraE_contributors,
			VEVMA_ia_reqm,
			VEVMA_ia_reqm_hbond,
			VEVMA_ia_atom1_types_reqm,
			VEVMA_ia_atom2_types_reqm,
			VEVMA_ia_VWpars_AC,
			VEVMA_ia_VWpars_BD,
			VEVMA_ia_dspars_S,
			VEVMA_ia_dspars_V,
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
