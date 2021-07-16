#include "auxiliary.h"
#include "rand_gen.h"
#include "calc_pc.h"
#include "energy_and_gradient.h"

// Stopping criterion from Solis-Wets
//#define ADADELTA_AUTOSTOP

// ADADELTA parameters
#define RHO 0.8f
#define EPSILON 1e-2f

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

void ls_ad(
			uint				DockConst_max_num_of_iters,
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
    printf("LS: DockConst_max_num_of_iters: %u\n",		DockConst_max_num_of_iters);
    printf("LS: DockConst_num_of_genes: %u\n",		    DockConst_num_of_genes);
    printf("LS: pop_size: %u\n",		                pop_size);
    printf("LS: DockConst_rotbondlist_length: %u\n",    DockConst_rotbondlist_length);
    printf("LS: Host_RunId: %u\n",                      Host_RunId);
    printf("LS: DockConst_smooth: %f\n",                DockConst_smooth);
    printf("LS: DockConst_num_of_intraE_contributors: %u\n",    DockConst_num_of_intraE_contributors);
    printf("LS: DockConst_grid_spacing: %f\n",          DockConst_grid_spacing);
    printf("LS: DockConst_num_of_atypes: %u\n",         DockConst_num_of_atypes);
    printf("LS: DockConst_coeff_elec: %f\n",            DockConst_coeff_elec);
    printf("LS: DockConst_qasp: %f\n",                  DockConst_qasp);
    printf("LS: DockConst_coeff_desolv: %f\n",          DockConst_coeff_desolv);
    printf("LS: DockConst_xsz: %u\n",                   DockConst_xsz);
    printf("LS: DockConst_ysz: %u\n",                   DockConst_ysz);
    printf("LS: DockConst_zsz: %u\n",                   DockConst_zsz);
    printf("LS: DockConst_num_of_atoms: %u\n",          DockConst_num_of_atoms);
    printf("LS: DockConst_gridsize_x_minus1: %f\n",     DockConst_gridsize_x_minus1);
    printf("LS: DockConst_gridsize_y_minus1: %f\n",     DockConst_gridsize_y_minus1);
    printf("LS: DockConst_gridsize_z_minus1: %f\n",     DockConst_gridsize_z_minus1);
    printf("LS: Host_mul_tmp2: %u\n",                   Host_mul_tmp2);
    printf("LS: Host_mul_tmp3: %u\n",                   Host_mul_tmp3);
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

    // Update vector, i.e., "delta"
    // It is added to the genotype to create the next genotype.
    // E.g., in steepest descent,: delta = -1.0 * stepsize * gradient
    float delta[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

    // Squared updates E[dx^2]
    float square_delta[ACTUAL_GENOTYPE_LENGTH][MAX_POPSIZE];

    // Initializing vectors
    for (uint i = 0; i < DockConst_num_of_genes; i++) {
        for (uint j = 0; j < pop_size; j++) {
            gradient[i][j]          = 0.0f;
            square_gradient[i][j]   = 0.0f;
            delta[i][j]             = 0.0f;
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

#ifdef ADADELTA_AUTOSTOP
	float rho = 1.0f;
	uint cons_succ = 0;
	uint cons_fail = 0;
#endif

    // TODO: add other blocks for ADADELTA autostop

    // Perfoming ADADELTA iterations

	// Iteration counter for the ADADELTA minimizer
	uint iteration_cnt = 0;

	// --------------------------------------------------------------------------
    // The termination criteria is based on a maximum number of iterations,
    // and the minimum step size allowed for single-precision FP numbers
    // (IEEE-754 single float has a precision of about 6 decimal digits)
    do {

#ifdef PRINT_ALL_LS_AD
	printf("LS_ADADELTA: iteration_cnt: %u\n", iteration-cnt);
#endif
		// TODO
		// Calculating energy and gradients
        float energy_ie[MAX_POPSIZE];
        float energy_ia[MAX_POPSIZE];

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
			pop_size /*active_pop_size*/,    // TODO: FIXME
			genotype /*entity_possible_new_genotype*/,   // TODO: FIXME
			local_coords_x,
			local_coords_y,
			local_coords_z
		);
        energy_and_gradient(
            genotype,
            energy_ie,
            energy_ia,
            local_coords_x,
            local_coords_y,
            local_coords_z,
            gradient,
            DockConst_num_of_genes,
            pop_size,
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

        // TODO: fix usage of j
        for (uint j = 0; j < pop_size; j++) {
            energy[j] = energy_ia[j] + energy_ie[j];

#ifdef PRINT_ALL_LS_AD
            printf("energy[%i]: %f\n", j, energy[j]);
#endif

        }

		for (uint i = 0; i < DockConst_num_of_genes; i++) {

            // TODO: fix usage of j
            for (uint j = 0; j < pop_size; j++) {
                if (energy[j] < best_energy[j]) {
                    best_genotype[i][j] = genotype[i][j];
                }

                // Accummulating gradient^2 (Eq.8 in paper)
                // square_gradient corresponds to E[g^2]
                square_gradient[i][j] = RHO * square_gradient[i][j] + (1.0f - RHO) * gradient[i][j] * gradient[i][j];

                // Computing update (Eq.9 in paper)
                float tmp_div = (square_delta[i][j] + EPSILON) / (square_gradient[i][j] + EPSILON);
                delta[i][j] = -1.0f * gradient[i][j] * esa_sqrt(tmp_div);

                // Accummulating update^2
                // square_delta corresponds to E[dx^2]
                square_delta[i][j] = RHO * square_delta[i][j] + (1.0f - RHO) * delta[i][j] * delta[i][j];

                // Applying update
                genotype[i][j] = genotype[i][j] + delta[i][j];
            } // End j Loop (over individuals)
		}

#ifdef PRINT_ALL_LS_AD
	    printf("\n%s\n", "----------------------------------------------------------");
	    printf("%13s %20s %15s %15s %15s\n", "gene", "sq_grad", "delta", "sq_delta", "new.genotype");
	    for (uint i = 0; i < DockConst_num_of_genes; i++) {
		    printf("%13u %20.6f %15.6f %15.6f %15.6f\n", i, square_gradient[i], delta[i], square_delta[i], genotype[i]);
        }
#endif

        // TODO: fix usage of j
        for (uint j = 0; j < pop_size; j++) {
            // Updating number of ADADELTA iterations (energy evaluations)
            if (energy[j] < best_energy[j]) {
                best_energy[j] = energy[j];

        #ifdef ADADELTA_AUTOSTOP
                cons_succ++;
                cons_fail = 0;
        #endif
            }
        #ifdef ADADELTA_AUTOSTOP
            else {
                cons_succ = 0;
                cons_fail++;
            }
        #endif

            iteration_cnt = iteration_cnt + 1;

        #ifdef ADADELTA_AUTOSTOP
            if (cons_succ >=4) {
                rho *= LS_EXP_FACTOR;
                cons_succ = 0;
            }
            else {
                if (cons_fail >= 4) {
                    rho *= LS_CONT_FACTOR;
                    cons_fail = 0;
                }
            }
        #endif
        } // End j Loop (over individuals)

#ifdef ADADELTA_AUTOSTOP
	} while ((iteration_cnt < DockConst_max_num_of_iters) && (rho > 0.01f));
#else
    } while (iteration_cnt < DockConst_max_num_of_iters);
#endif
	// --------------------------------------------------------------------------

	// TODO: double check if ADGPU uses two types of mapping functions
	// Mapping angles
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
        // TODO: fix usage of j
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
        // TODO: fix usage of j
        for (uint j = 0; j < pop_size; j++) {
		    in_out_genotype[i][j] = best_genotype[i][j];
        } // End j Loop (over individuals)
	}

	// Updating energy
    for (uint j = 0; j < pop_size; j++) {
	    in_out_energy[j] = best_energy[j];
    } // End j Loop (over individuals)

	// Updating evals
	*out_eval = iteration_cnt;

#ifdef PRINT_ALL_LS_AD
	printf("\n");
	printf("Finishing <local search>: <ADADELTA>\n");
	printf("\n");
#endif
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
