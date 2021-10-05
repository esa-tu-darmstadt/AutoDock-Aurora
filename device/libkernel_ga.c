/*

AutoDock-Aurora, a vectorized implementation of AutoDock 4.2 for the NEC SX-Aurora TSUBASA
Copyright (C) 2021 TU Darmstadt, Embedded Systems and Applications Group, Germany. All rights reserved.
For some of the code, Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

AutoDock is a Trade Mark of the Scripps Research Institute.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include "auxiliary.h"
#include "rand_gen.h"
#include "device_args.h"
#include "lga.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

// --------------------------------------------------------------------------
// Lamarckian Genetic-Algorithm (GA): GA + LS (Local Search)
// --------------------------------------------------------------------------
uint64_t libkernel_ga(
                       struct device_args *da
                     )
{
	uint DockConst_pop_size = da->DockConst_pop_size;
	uint Host_num_of_runs = da->Host_num_of_runs;
	int ve_proc_id = da->ve_proc_id;
	int ve_num_procs = da->ve_num_procs;

	// Initializing random generator with seed passed from host
	randf_vec_init(&(da->dockpars_prng_states[ve_proc_id * DockConst_pop_size]), DockConst_pop_size);

	int ve_num_of_runs = (Host_num_of_runs + ve_num_procs - 1) / ve_num_procs;
	uint run_cnt_start = ve_proc_id * ve_num_of_runs;
	uint run_cnt_end = MIN((ve_proc_id + 1) * ve_num_of_runs, Host_num_of_runs);

#pragma omp parallel for schedule(static, 1)
	for (uint run_cnt = run_cnt_start; run_cnt < run_cnt_end; run_cnt++) {

#ifdef PRINT_ALL_KRNL
		printf(" %u", run_cnt+1); 
		fflush(stdout);
#endif

		// Values changing every LGA run
		uint uint_run_cnt  = run_cnt; // - run_cnt_start;
		uint Host_Offset_Pop = uint_run_cnt * DockConst_pop_size * ACTUAL_GENOTYPE_LENGTH;
		uint Host_Offset_Ene = uint_run_cnt * DockConst_pop_size;

		lga(
			da->PopulationCurrentInitial,
			da->PopulationCurrentFinal,
			da->EnergyCurrent,
			da->Evals_performed,
			da->Gens_performed,
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
			da->DockConst_num_of_genes,
			// PC
			da->PC_rotlist,
			da->PC_ref_coords_x, // TODO: merge them into a single one?
			da->PC_ref_coords_y,
			da->PC_ref_coords_z,
			da->PC_rotbonds_moving_vectors,
			da->PC_rotbonds_unit_vectors,
			da->PC_ref_orientation_quats,
			da->DockConst_rotbondlist_length,
			// IA
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
			// IE
			da->Fgrids,
			da->DockConst_xsz,
			da->DockConst_ysz,
			da->DockConst_zsz,
			da->DockConst_num_of_atoms,
			da->DockConst_gridsize_x_minus1,
			da->DockConst_gridsize_y_minus1,
			da->DockConst_gridsize_z_minus1,
			da->Host_mul_tmp2,
			da->Host_mul_tmp3,
			// LS
			da->lsmet,
			// LS-SW
			da->DockConst_max_num_of_iters,
			da->DockConst_rho_lower_bound,
			da->DockConst_base_dmov_mul_sqrt3,
			da->DockConst_base_dang_mul_sqrt3,
			da->DockConst_cons_limit,
			// LS-AD
			da->GRAD_rotbonds,
			da->GRAD_rotbonds_atoms,
			da->GRAD_num_rotating_atoms_per_rotbond,
			da->GRAD_angle,
			da->GRAD_dependence_on_theta,
			da->GRAD_dependence_on_rotangle,
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
