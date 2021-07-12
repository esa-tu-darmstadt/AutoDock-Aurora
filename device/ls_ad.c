#include "auxiliary.h"
#include "rand_gen.h"
#include "calc_pc.h"
#include "energy_and_gradient.h"

// Stopping criterion from Solis-Wets
//#define ADADELTA_AUTOSTOP

// ADADELTA parameters
#define RHO 0.8f
#define EPSILON 1e-2 // TODO: make sure it is in SP FP

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

void ls_ad(
			uint				DockConst_max_num_of_iters,
            uchar               DockConst_num_of_genes, // ADGPU defines it as int
            uint                pop_size,
            float               in_out_genotype[][MAX_POPSIZE], // ADGPU: __global float* restrict dockpars_conformations_next,
            float*  restrict    in_out_energy, // ADGPU: __global float* restrict dockpars_energies_next,
            uint*   restrict    out_eval, // ADGPU: __global int* restrict dockpars_evals_of_new_entities,

    //uint DockConst_gridsize_x_times_y, // FIXME: present in ADPGU, but not in SW
    //uint DockConst_gridsize_x_times_y_times_z, // FIXME: present in ADGPU, but not in SW

    // pc
    const   int*    restrict    PC_rotlist,
    const   float*  restrict    PC_ref_coords_x,
    const   float*  restrict    PC_ref_coords_y,
    const   float*  restrict    PC_ref_coords_z,
    const   float*  restrict    PC_rotbonds_moving_vectors,
    const   float*  restrict    PC_rotbonds_unit_vectors,
    const   float*  restrict    PC_ref_orientation_quats,
            uint                DockConst_rotbondlist_length,
            uint                Host_RunId,
    // ia
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
    // ie
    const   float*  restrict    IE_Fgrids,
            uchar               DockConst_xsz,
            uchar               DockConst_ysz,
            uchar               DockConst_zsz,
            uint                DockConst_num_of_atoms,    // TODO FIXME: SW defines it as uchar
            float               DockConst_gridsize_x_minus1, // TODO FIXME: ADGPU defines it as int
            float               DockConst_gridsize_y_minus1,
            float               DockConst_gridsize_z_minus1,
            uint                Host_mul_tmp2,
            uint                Host_mul_tmp3,
    // gradients
    const   int*    restrict    GRAD_rotbonds,
    const   int*    restrict    GRAD_rotbonds_atoms,
    const   int*    restrict    GRAD_num_rotating_atoms_per_rotbond,

    const   float*  restrict    GRAD_angle,
    const   float*  restrict    GRAD_dependence_on_theta,
    const   float*  restrict    GRAD_dependence_on_rotangle
) {
	// Genotype and its energy
	float genotype[ACTUAL_GENOTYPE_LENGTH];
	float energy;

    // Partial results of the gradient step
	float gradient[ACTUAL_GENOTYPE_LENGTH];

	// Energy may go up, so we keep track of the best energy ever calculated.
	// Then, we return the genotype corresponding to the best
	// observed energy, i.e., "best genotype"
    float best_energy;
    float best_genotype[ACTUAL_GENOTYPE_LENGTH];

    // Gradient of the intermolecular energy per each ligand atom
    // Also used to store the accummulated gradient per each ligand atom
    float gradient_inter_x[MAX_NUM_OF_ATOMS];
    float gradient_inter_y[MAX_NUM_OF_ATOMS];
    float gradient_inter_z[MAX_NUM_OF_ATOMS];

    // Gradient of the intramolecular energy per each ligand atom
    float gradient_intra_x[MAX_NUM_OF_ATOMS];
    float gradient_intra_y[MAX_NUM_OF_ATOMS];
    float gradient_intra_z[MAX_NUM_OF_ATOMS];

    // Ligand-atom position and partial energies
    float local_coords_x[MAX_NUM_OF_ATOMS]; //__local float4 calc_coords[MAX_NUM_OF_ATOMS];
    float local_coords_y[MAX_NUM_OF_ATOMS];
    float local_coords_z[MAX_NUM_OF_ATOMS];
    float local_coords_w[MAX_NUM_OF_ATOMS];

/*
	float local_coords_x[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_y[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
	float local_coords_z[MAX_NUM_OF_ATOMS][MAX_POPSIZE];
*/
	for (uint i = 0; i < DockConst_num_of_atoms; i++) {
/*
		for (uint j = 0; j < pop_size; j++) {
			local_coords_x[i][j] = 0.0f;
			local_coords_y[i][j] = 0.0f;
			local_coords_z[i][j] = 0.0f;
		}
*/
        local_coords_x[i] = 0.0f;
        local_coords_y[i] = 0.0f;
        local_coords_z[i] = 0.0f;

        gradient_inter_x[i] = 0.0f;
        gradient_inter_y[i] = 0.0f;
        gradient_inter_z[i] = 0.0f;

        gradient_intra_x[i] = 0.0f;
        gradient_intra_y[i] = 0.0f;
        gradient_intra_z[i] = 0.0f;
	}

    //float partial_energies[/*NUM_OF_THREADS_PER_BLOCK*/]; // TODO: remove it?

    // Squared gradients E[g^2]
    float square_gradient[ACTUAL_GENOTYPE_LENGTH];

    // Update vector, i.e., "delta"
    // It is added to the genotype to create the next genotype.
    // E.g., in steepest descent,: delta = -1.0 * stepsize * gradient
    float delta[ACTUAL_GENOTYPE_LENGTH];

    // Squared updates E[dx^2]
    float square_delta[ACTUAL_GENOTYPE_LENGTH];

    // Initializing vectors
    for (uint i= 0; i < DockConst_num_of_genes; i++) {
		gradient[i] = 0.0f;
        square_gradient[i] = 0.0f;
        delta[i] = 0.0f;
        square_delta[i] = 0.0f;
		genotype[i] = in_out_genotype[i]; // TODO: add 2dimension
        best_genotype[i] = in_out_genotype[i]; // TODO: add 2dimension
    }

	// Initializing energy
	energy = *in_out_energy; // TODO: make sure this read is correct

    // Initializing best_energy
    best_energy = INFINITY; // TODO: check correctness

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

        float energy_ie = 0; // TODO: will use eventually: float energy_ie[MAX_POPSIZE];
        float energy_ia = 0; // TODO: will use eventually: float energy_ia[MAX_POPSIZE];

		energy_and_gradient(
            genotype,           // TODO: verify it is correct!
            &energy_ie,
            &energy_ia,
            local_coords_x,
            local_coords_y,
            local_coords_z,
            gradient_inter_x,
            gradient_inter_y,
            gradient_inter_z,
            gradient_intra_x,
            gradient_intra_y,
            gradient_intra_z,
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

		for (uint i = 0; i < DockConst_num_of_genes; i++) {

			if (energy < best_energy) {
				best_genotype[i] = genotype[i];
			}

			// Accummulating gradient^2 (Eq.8 in paper)
			// square_gradient corresponds to E[g^2]
			square_gradient[i] = RHO * square_gradient[i] + (1.0f - RHO) * gradient[i] * gradient[i];

            // Computing update (Eq.9 in paper)
            float tmp_div = (square_delta[i] + EPSILON) / (square_gradient[i] + EPSILON);
            delta[i] = -1.0f * gradient[i] * esa_sqrt(tmp_div);

            // Accummulating update^2
            // square_delta corresponds to E[dx^2]
            square_delta[i] = RHO * square_delta[i] + (1.0f - RHO) * delta[i] * delta[i];

            // Applying update
            genotype[i] = genotype[i] + delta[i];

		}

#ifdef PRINT_ALL_LS_AD
	printf("\n%s\n", "----------------------------------------------------------");
	printf("%13s %20s %15s %15s %15s\n", "gene", "sq_grad", "delta", "sq_delta", "new.genotype");
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
		printf("%13u %20.6f %15.6f %15.6f %15.6f\n", i, square_gradient[i], delta[i], square_delta[i], genotype[i]);
	}
#endif

	// Updating number of ADADELTA iterations (energy evaluations)
	if (energy < best_energy) {
		best_energy = energy;

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

#ifdef ADADELTA_AUTOSTOP
	} while ((iteration_cnt < DockConst_max_num_of_iters) && (rho > 0.01f));
#else
    } while (iteration_cnt < DockConst_max_num_of_iters);
#endif
	// --------------------------------------------------------------------------

	// TODO: double check if ADGPU uses two types of mapping functions
	// Mapping angles
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
		if (i > 2) {
			if (i == 4) {
				best_genotype[i] = map_angle_180(best_genotype[i]);
			} else {
				best_genotype[i] = map_angle_360(best_genotype[i]);
			}
		}
	}

	// TODO: vectorize it
	// Updating old offspring in population
	for (uint i = 0; i < DockConst_num_of_genes; i++) {
		in_out_genotype[i] = best_genotype[i];
	}

	// TODO: vectorize it
	// Updating energy
	in_out_energy = best_energy;

	// Updating evals
	*out_eval = iteration_cnt;
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
