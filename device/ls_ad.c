#include "auxiliary.h"
#include "rand_gen.h"
#include "calc_pc.h"
#include "energy_and_gradient.h"

// Stopping criterion from Solis-Wets
//#define ADADELTA_AUTOSTOP

// ADADELTA parameters
#define RHO 0.8f
#define EPSILON 1e-2f

// Enable for debugging ADADELTA using a defined initial genotype
//#define DEBUG_ADADELTA_INITIAL_2BRT

#ifdef DEBUG_ADADELTA_INITIAL_2BRT
#include "energy_ia.h"
#include "energy_ie.h"
#endif

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

void ls_ad(
			ushort				DockConst_max_num_of_iters,
            uchar               DockConst_num_of_genes, // TODO: ADGPU defines it as int
            uint                pop_size,
            float               in_out_genotype[][MAX_POPSIZE], // TODO: ADGPU: __global float* restrict dockpars_conformations_next,
            float*  restrict    in_out_energy, // TODO: ADGPU: __global float* restrict dockpars_energies_next,
            uint*   restrict    out_eval, // TODO: ADGPU: __global int* restrict dockpars_evals_of_new_entities,

    //uint DockConst_gridsize_x_times_y, // TODO: present in ADPGU, but not in SW
    //uint DockConst_gridsize_x_times_y_times_z, // TODO: present in ADGPU, but not in SW

    // PC
    const   int*    restrict    PC_rotlist,
    const   float*  restrict    PC_ref_coords_x,
    const   float*  restrict    PC_ref_coords_y,
    const   float*  restrict    PC_ref_coords_z,
    const   float*  restrict    PC_rotbonds_moving_vectors,
    const   float*  restrict    PC_rotbonds_unit_vectors,
    const   float*  restrict    PC_ref_orientation_quats,
            uint                DockConst_rotbondlist_length,
            uint                Host_RunId,
    // IA
    const   float*  restrict    IA_IE_atom_charges,
    const   int*    restrict    IA_IE_atom_types,
    const   int*    restrict    IA_intraE_contributors,
    const   float*  restrict    IA_reqm,
    const   float*  restrict    IA_reqm_hbond,
    const   uint*   restrict    IA_atom1_types_reqm,
    const   uint*   restrict    IA_atom2_types_reqm,
    const   float*  restrict    IA_VWpars_AC,
    const   float*  restrict    IA_VWpars_BD,
    const   float*  restrict    IA_dspars_S,
    const   float*  restrict    IA_dspars_V,
            float               DockConst_smooth,
            uint                DockConst_num_of_intraE_contributors,
            float               DockConst_grid_spacing,
            uint                DockConst_num_of_atypes,
            float               DockConst_coeff_elec,
            float               DockConst_qasp,
            float               DockConst_coeff_desolv,
    // IE
    const   float*  restrict    IE_Fgrids,
            uchar               DockConst_xsz,
            uchar               DockConst_ysz,
            uchar               DockConst_zsz,
            uchar               DockConst_num_of_atoms,
            float               DockConst_gridsize_x_minus1, // TODO FIXME: ADGPU defines it as int
            float               DockConst_gridsize_y_minus1,
            float               DockConst_gridsize_z_minus1,
            uint                Host_mul_tmp2,
            uint                Host_mul_tmp3,
    // Gradients
    const   int*    restrict    GRAD_rotbonds,
    const   int*    restrict    GRAD_rotbonds_atoms,
    const   int*    restrict    GRAD_num_rotating_atoms_per_rotbond,
    const   float*  restrict    GRAD_angle,
    const   float*  restrict    GRAD_dependence_on_theta,
    const   float*  restrict    GRAD_dependence_on_rotangle
) {
#ifdef PRINT_ALL_LS_AD
    printf("\n");
    printf("Starting <local search>: <ADADELTA> ... \n");
    printf("\n");
    printf("%-50s %u\n", "AD / DockConst_max_num_of_iters: ",      DockConst_max_num_of_iters);
    printf("%-50s %u\n", "AD / DockConst_num_of_genes: ",          DockConst_num_of_genes);
    printf("%-50s %u\n", "AD / pop_size: ",		                   pop_size);
    printf("%-50s %u\n", "AD / DockConst_rotbondlist_length: ",    DockConst_rotbondlist_length);
    printf("%-50s %u\n", "AD / Host_RunId: ",                      Host_RunId);
    printf("%-50s %f\n", "AD / DockConst_smooth: ",                DockConst_smooth);
    printf("%-50s %u\n", "AD / DockConst_num_of_intraE_contributors: ",    DockConst_num_of_intraE_contributors);
    printf("%-50s %f\n", "AD / DockConst_grid_spacing: ",          DockConst_grid_spacing);
    printf("%-50s %u\n", "AD / DockConst_num_of_atypes: ",         DockConst_num_of_atypes);
    printf("%-50s %f\n", "AD / DockConst_coeff_elec: ",            DockConst_coeff_elec);
    printf("%-50s %f\n", "AD / DockConst_qasp: ",                  DockConst_qasp);
    printf("%-50s %f\n", "AD / DockConst_coeff_desolv: ",          DockConst_coeff_desolv);
    printf("%-50s %u\n", "AD / DockConst_xsz: ",                   DockConst_xsz);
    printf("%-50s %u\n", "AD / DockConst_ysz: ",                   DockConst_ysz);
    printf("%-50s %u\n", "AD / DockConst_zsz: ",                   DockConst_zsz);
    printf("%-50s %u\n", "AD / DockConst_num_of_atoms: ",          DockConst_num_of_atoms);
    printf("%-50s %f\n", "AD / DockConst_gridsize_x_minus1: ",     DockConst_gridsize_x_minus1);
    printf("%-50s %f\n", "AD / DockConst_gridsize_y_minus1: ",     DockConst_gridsize_y_minus1);
    printf("%-50s %f\n", "AD / DockConst_gridsize_z_minus1: ",     DockConst_gridsize_z_minus1);
    printf("%-50s %u\n", "AD / Host_mul_tmp2: ",                   Host_mul_tmp2);
    printf("%-50s %u\n", "AD / Host_mul_tmp3: ",                   Host_mul_tmp3);
#endif

	// Genotype and its energy
	float genotype[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];
	float energy[MAX_POPSIZE];

    // Partial results of the gradient step
	float gradient[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

	// Energy may go up, so we keep track of the best energy ever calculated.
	// Then, we return the genotype corresponding to the best
	// observed energy, i.e., "best genotype"
    float best_energy[MAX_POPSIZE];
    float best_genotype[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

    // Ligand-atom position and partial energies
	float local_coords_x[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_y[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_z[MAX_NUM_OF_ATOMS][MAX_POPSIZE];

	for (uint i = 0; i < DockConst_num_of_atoms; i++) {
		for (uint j = 0; j < pop_size; j++) {
			local_coords_x[i][j] = 0.0f;
			local_coords_y[i][j] = 0.0f;
			local_coords_z[i][j] = 0.0f;
		}
	}

    // Squared gradients E[g^2]
    float square_gradient[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

    // Transformed into scalar
    // float delta[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

    // TODO: convert to scalar
    // Squared updates E[dx^2]
    float square_delta[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

    // Initializing vectors
    for (uint i = 0; i < DockConst_num_of_genes; i++) {
        for (uint j = 0; j < pop_size; j++) {
            // gradient[i][j]          = 0.0f; // Initialized in <energy_and_gradient()>
            square_gradient[i][j]   = 0.0f;
            square_delta[i][j]      = 0.0f;
            genotype[i][j]          = in_out_genotype[i][j];
            best_genotype[i][j]     = in_out_genotype[i][j];
        } // End j Loop (over individuals)
    }

	// Initializing energy & best_energy
    for (uint j = 0; j < pop_size; j++) {
        energy[j] = in_out_energy[j];
        best_energy[j] = INFINITY;
    } // End j Loop (over individuals)

    float energy_ie_ad[MAX_POPSIZE];
    float energy_ia_ad[MAX_POPSIZE];

    // WARNING: hardcoded definitions overwrite assignments from LGA.
    // This means that when <DEBUG_ADADELTA_INITIAL_2BRT> is defined,
	// then, the LGA genotype is not used (only for debugging purposes!)
#ifdef DEBUG_ADADELTA_INITIAL_2BRT
    for (uint i = 0; i < DockConst_num_of_genes; i++) {
        for (uint j = 0; j < pop_size; j++) {
            // 2brt
            genotype[0][j]  = 24.093334f;
            genotype[1][j]  = 24.658667f;
            genotype[2][j]  = 24.210667f;
            genotype[3][j]  = 50.0f;
            genotype[4][j]  = 50.0f;
            genotype[5][j]  = 50.0f;
            genotype[6][j]  = 0.0f;
            genotype[7][j]  = 0.0f;
            genotype[8][j]  = 0.0f;
            genotype[9][j]  = 0.0f;
            genotype[10][j] = 0.0f;
            genotype[11][j] = 0.0f;
            genotype[12][j] = 0.0f;
            genotype[13][j] = 0.0f;
            genotype[14][j] = 0.0f;
            genotype[15][j] = 0.0f;
            genotype[16][j] = 0.0f;
            genotype[17][j] = 0.0f;
            genotype[18][j] = 0.0f;
            genotype[19][j] = 0.0f;
            genotype[20][j] = 0.0f;
        }
    }

    calc_pc(
		PC_rotlist,
        PC_ref_coords_x,
        PC_ref_coords_y,
        PC_ref_coords_z,
        PC_rotbonds_moving_vectors,
        PC_rotbonds_unit_vectors,
        PC_ref_orientation_quats,
        DockConst_rotbondlist_length,
        DockConst_num_of_genes,
        Host_RunId,
        pop_size,
        genotype,
        local_coords_x,
        local_coords_y,
        local_coords_z
	);
	energy_ia(
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
        pop_size,
        energy_ia_ad,
        local_coords_x,
        local_coords_y,
        local_coords_z
	);
	energy_ie(
        IE_Fgrids,
        IA_IE_atom_charges,
        IA_IE_atom_types,
        DockConst_xsz,
        DockConst_ysz,
        DockConst_zsz,
        DockConst_num_of_atoms,
        DockConst_gridsize_x_minus1,
        DockConst_gridsize_y_minus1,
        DockConst_gridsize_z_minus1,
        Host_mul_tmp2,
        Host_mul_tmp3,
        pop_size,
        energy_ie_ad,
        local_coords_x,
        local_coords_y,
        local_coords_z
	);

    for (uint j = 0; j < pop_size; j++) {
		energy[j] = energy_ia_ad[j] + energy_ie_ad[j];
	}

	printf("\n");
	printf("%20s %u \n", "hardcoded genotype: ", 0);
	printf("%20s %.6f\n", "initial energy: ", energy[0]);
#endif

#ifdef ADADELTA_AUTOSTOP
	float rho = 1.0f;
	uint cons_succ = 0;
	uint cons_fail = 0;
#endif

    // TODO: add other blocks for ADADELTA autostop

    // Perfoming ADADELTA iterations

	// Iteration counter for the ADADELTA minimizer
	//uint iteration_cnt = 0;
    uint iteration_cnt[MAX_POPSIZE];
    uint iteration_compr[MAX_POPSIZE];

    uint active_pop_size; // Counts compressed list of genotypes
    uint active_idx[MAX_POPSIZE]; // Index that corresponds to slot in compressed list

    uint LS_eval = 0;

    uint ls_is_active[MAX_POPSIZE]; // Filter for individuals that are still being iterated
    uint num_active_ls = pop_size;

    for (uint j = 0; j < pop_size; j++) {
        iteration_cnt[j] = 0;
        ls_is_active[j] = 1;
    }

	// --------------------------------------------------------------------------
    // The termination criteria is based on a maximum number of iterations,
    // and the minimum step size allowed for single-precision FP numbers
    // (IEEE-754 single float has a precision of about 6 decimal digits)
    do {
        // Compressed list of active indexes
        active_pop_size = 0;
        for (uint j = 0; j < pop_size; j++) {
            if (ls_is_active[j]) {
                active_idx[active_pop_size] = j;
                iteration_compr[active_pop_size] = iteration_cnt[j];
                active_pop_size++;
            }
        }

        for (uint jj = 0; jj < active_pop_size; jj++) {
            iteration_compr[jj]++;
        }

        // Calculating energy and gradients
        calc_pc (
            PC_rotlist,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			DockConst_rotbondlist_length,
			DockConst_num_of_genes,
			Host_RunId,
			active_pop_size,
			genotype /*entity_possible_new_genotype*/,   // TODO: FIXME
			local_coords_x,
			local_coords_y,
			local_coords_z
		);
        energy_and_gradient(
            genotype,
            energy_ie_ad,
            energy_ia_ad,
            local_coords_x,
            local_coords_y,
            local_coords_z,
            gradient,
            DockConst_num_of_genes,
            active_pop_size,
            //ie
            IE_Fgrids,
            IA_IE_atom_charges,
            IA_IE_atom_types,
            DockConst_xsz,
            DockConst_ysz,
            DockConst_zsz,
            DockConst_num_of_atoms,
            DockConst_gridsize_x_minus1,
            DockConst_gridsize_y_minus1,
            DockConst_gridsize_z_minus1,
            Host_mul_tmp2,
            Host_mul_tmp3,
            //ia
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
            //gradients
            GRAD_rotbonds,
            GRAD_rotbonds_atoms,
            GRAD_num_rotating_atoms_per_rotbond,
            GRAD_angle,
            GRAD_dependence_on_theta,
            GRAD_dependence_on_rotangle
        );

        for (uint jj = 0; jj < active_pop_size; jj++) {
            energy[jj] = energy_ia_ad[jj] + energy_ie_ad[jj];

#ifdef PRINT_ALL_LS_AD
            if (jj == 0) {
                printf("\n");
                printf("%-10s %-10.6f \n", "intra: ",  energy_ia_ad[jj]);
                printf("%-10s %-10.6f \n", "grids: ",  energy_ie_ad[jj]);
                printf("%-10s %-10.6f \n", "Energy: ", energy[jj]);
                //printf("\t energy[%i]: %.3f\n", jj, energy[jj]);

                for(uint i = 0; i < DockConst_num_of_genes; i++) {
                    if (i == 0) {
                        printf("\n%s\n", "----------------------------------------------------------");
                        printf("%13s %13s %5s %15s %15s\n", "gene_id", "gene.value", "|", "gene.grad", "(autodockdevpy units)");
                    }
                    printf("%13u %13.6f %5s %15.6f %15.6f\n", i, genotype[i][jj], "|", gradient[i][jj], (i<3)? (gradient[i][jj]/0.375f):(gradient[i][jj]*180.0f/PI_FLOAT));
                }

                for(uint i = 0; i < DockConst_num_of_atoms; i++) {
                    if (i == 0) {
                        printf("\n%s\n", "----------------------------------------------------------");
                        printf("%s\n", "Coordinates");
                        printf("%12s %12s %12s %12s\n", "atom_id", "coords.x", "coords.y", "coords.z");
                    }
                    printf("%12u %12.6f %12.6f %12.6f\n", i, local_coords_x[i][jj], local_coords_y[i][jj], local_coords_z[i][jj]);
                }
                printf("\n");

            }
#endif
        } // End jj Loop (over active individuals)

        // Updating LS energy-evaluation count
		LS_eval += active_pop_size;

        uint energy_lower[MAX_POPSIZE];
        for (uint jj = 0; jj < active_pop_size; jj++) {
            uint j = active_idx[jj];
            if (energy[jj] < best_energy[j]) {
                energy_lower[jj] = 1;
            } else {
                energy_lower[jj] = 0;
            }
        } // End jj Loop (over active individuals)

        for (uint i = 0; i < DockConst_num_of_genes; i++) {
            for (uint jj = 0; jj < active_pop_size; jj++) {
                uint j = active_idx[jj];
                if (energy_lower[jj]) {
                    best_genotype[i][j] = genotype[i][jj];
                }
            } // End jj Loop (over active individuals)
        }

        for (uint i = 0; i < DockConst_num_of_genes; i++) {
            for (uint jj = 0; jj < active_pop_size; jj++) {
                uint j = active_idx[jj];

                // Accummulating gradient^2 (Eq.8 in paper)
                // square_gradient corresponds to E[g^2]
                square_gradient[i][j] = RHO * square_gradient[i][j] + (1.0f - RHO) * gradient[i][j] * gradient[i][j];

                // Computing update (Eq.9 in paper)
                float tmp_div = (square_delta[i][j] + EPSILON) / (square_gradient[i][j] + EPSILON);

                // Update or "delta"
                // It is added to the genotype to create the next genotype.
                // E.g., in steepest descent,: delta = -1.0 * stepsize * gradient
                float delta = -1.0f * gradient[i][j] * /*esa_sqrt*/sqrtf(tmp_div);

                // Accummulating update^2
                // square_delta corresponds to E[dx^2]
                square_delta[i][j] = RHO * square_delta[i][j] + (1.0f - RHO) * delta * delta;

                // Applying update
                genotype[i][jj] = genotype[i][jj] + delta;
            } // End jj Loop (over active individuals)

#ifdef PRINT_ALL_LS_AD
            if (i == 0) {
                printf("\n%s\n", "----------------------------------------------------------");
                printf("%13s %20s %15s %15s %15s\n", "gene", "sq_grad", "delta", "sq_delta", "new.genotype");
            }
            printf("%13u %15.6f %15.6f %15.6f %15.6f\n", i, square_gradient[i][0], delta[i][0], square_delta[i][0], genotype[i][0]);
#endif
		} // End Loop (over DockConst_num_of_genes)

#ifdef PRINT_ALL_LS_AD
	    printf("LS_ADADELTA: iteration_cnt: %u\n", iteration_cnt);
#endif

        for (uint jj = 0; jj < active_pop_size; jj++) {
            uint j = active_idx[jj];

            // Updating number of ADADELTA iterations (energy evaluations)
            if (energy_lower[jj]) {

#ifdef PRINT_ALL_LS_AD
                //printf("\t Improved ind %i: energy: %.3f, best_energy: %.3f\n", j, energy[j], best_energy[j]);
                printf("\t %-10s %i \t %-10s %15.3f \t %-10s %15.3f \t %-10s\n", "Ind ", j, "E: ", energy[j], "BE: ", best_energy[j], "(Improved!)");
#endif
                best_energy[j] = energy[jj];

#ifdef ADADELTA_AUTOSTOP
                cons_succ++;
                cons_fail = 0;
#endif

            } // First end: if (energy comparison)

#ifdef ADADELTA_AUTOSTOP
            else {
                cons_succ = 0;
                cons_fail++;
            } // Second end: if (energy comparison)

            if (cons_succ >= 4) {
                rho *= LS_EXP_FACTOR;
                cons_succ = 0;
            } else if (cons_fail >= 4) {
                rho *= LS_CONT_FACTOR;
                cons_fail = 0;
            }
#endif
        } // End jj Loop (over active individuals)

        num_active_ls = active_pop_size;
        for (uint jj = 0; jj < active_pop_size; jj++) {
            uint j = active_idx[jj];
            if (iteration_compr[jj] > DockConst_max_num_of_iters) {
                ls_is_active[j] = 0;
                num_active_ls--;
            }
            iteration_cnt[j] = iteration_compr[jj];
        } // End jj Loop (over active individuals)

#ifdef ADADELTA_AUTOSTOP
	} while ((iteration_cnt < DockConst_max_num_of_iters) && (rho > 0.01f));
#else
    //} while (iteration_cnt < DockConst_max_num_of_iters);
    } while (num_active_ls > 0);
#endif
	// --------------------------------------------------------------------------

	// TODO: double check if ADGPU uses two types of mapping functions
	// Mapping angles
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
        for (uint j = 0; j < pop_size; j++) {
            if (i > 2) {
                if (i == 4) {
                    best_genotype[i][j] = map_angle_180(best_genotype[i][j]);
                } else {
                    best_genotype[i][j] = map_angle_360(best_genotype[i][j]);
                }
            }
        } // End j Loop (over individuals)
	}

	// Updating old offspring in population
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
        for (uint j = 0; j < pop_size; j++) {
		    in_out_genotype[i][j] = best_genotype[i][j];
        } // End j Loop (over individuals)
	}

	// Updating energy
    for (uint j = 0; j < pop_size; j++) {
	    in_out_energy[j] = best_energy[j];
    } // End j Loop (over individuals)

	// Updating evals
	*out_eval = LS_eval;

#ifdef PRINT_ALL_LS_AD
	printf("\n");
	printf("Finishing <local search>: <ADADELTA>\n");
	printf("\n");
#endif
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
