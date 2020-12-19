#include "auxiliary.h"
#include "rand_gen.h"
#include "device_args.h"
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
                       struct device_args *da
                       )
#if 0
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
#endif
{
	uint DockConst_pop_size = da->DockConst_pop_size;
	uint Host_num_of_runs = da->Host_num_of_runs;

	// initialize random generator with seed passed from host
	randf_vec_init(&(da->dockpars_prng_states)[0], DockConst_pop_size);


#pragma omp parallel for schedule(static, 1)
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
			da->PopulationCurrentInitial,
			da->PopulationCurrentFinal,
			da->EnergyCurrent,
			da->Evals_performed,
			da->Gens_performed,
			da->dockpars_prng_states,
			da->DockConst_pop_size,
			da->DockConst_num_of_energy_evals,
			da->DockConst_num_of_generations,
			da->DockConst_tournament_rate,
			da->DockConst_mutation_rate,
			da->DockConst_abs_max_dmov,
			da->DockConst_abs_max_dang,
			da->Host_two_absmaxdmov,
			da->Host_two_absmaxdang,
			da->DockConst_crossover_rate,
			da->DockConst_num_of_lsentities,
			da->DockConst_num_of_genes,
			// pc
			da->PC_rotlist,
			da->PC_ref_coords_x,// TODO: merge them into a single one?
			da->PC_ref_coords_y,
			da->PC_ref_coords_z,
			da->PC_rotbonds_moving_vectors,
			da->PC_rotbonds_unit_vectors,
			da->PC_ref_orientation_quats,
			da->DockConst_rotbondlist_length,
			// ia
			da->IA_IE_atom_charges,
			da->IA_IE_atom_types,
			da->IA_intraE_contributors,
			da->IA_reqm,
			da->IA_reqm_hbond,
			da->IA_atom1_types_reqm,
			da->IA_atom2_types_reqm,
			da->IA_VWpars_AC,
			da->IA_VWpars_BD,
			da->IA_dspars_S,
			da->IA_dspars_V,
			da->DockConst_smooth,
			da->DockConst_num_of_intraE_contributors,
			da->DockConst_grid_spacing,
			da->DockConst_num_of_atypes,
			da->DockConst_coeff_elec,
			da->DockConst_qasp,
			da->DockConst_coeff_desolv,
			// ie
			da->Fgrids,
			da->DockConst_g1,
			da->DockConst_g2,
			da->DockConst_g3,
			da->DockConst_num_of_atoms,
			da->DockConst_gridsize_x_minus1,
			da->DockConst_gridsize_y_minus1,
			da->DockConst_gridsize_z_minus1,
			da->Host_mul_tmp2,
			da->Host_mul_tmp3,
			// ls
			da->DockConst_max_num_of_iters,
			da->DockConst_rho_lower_bound,
			da->DockConst_base_dmov_mul_sqrt3,
			da->DockConst_base_dang_mul_sqrt3,
			da->DockConst_cons_limit,
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
