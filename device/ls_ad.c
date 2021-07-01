#include "auxiliary.h"
#include "rand_gen.h"
#include "energy_ia.h"
#include "energy_ie.h"
#include "calc_pc.h"

// Stopping criterion from Solis-Wets
//#define ADADELTA_AUTOSTOP

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
    const   int*    restrict    IA_IA_atom_types,
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
            float               DockConst_grd_spacing,
            uint                DockConst_num_of_atypes,
            float               DockConst_coeff_elec,
            float               DockConst_qasp,
            float               DockConst_coeff_desolv,
    // ie
    const   float*  restrict    IE_Fgrids,
            uchar               DockConst_xsz,
            uchar               DockConst_ysz,
            uchar               DockConst_zsz,
            uint                DockConst_num_of_atoms,    // FIXME: SW defines it as uchar
            float               DockConst_gridsize_x_minus1, // FIXME: ADGPU defines it as int
            float               DockConst_gridsize_y_minus1,
            float               DockConst_gridsize_z_minus1,
            uint                Host_mul_tmp2,
            uint                Host_mul_tmp3
    // gradients
    const   int*    restrict    GRAD_rotbonds,
    const   int*    restrict    GRAD_rotbonds_atoms,
    const   int*    restrict    GRAD_num_rotating_atoms_per_rotbond,

    const   float*  restrict    GRAD_angle,
    const   float*  restrict    GRAD_dependence_on_theta,
    const   float*  restrict    GRAD_dependence_on_rotangle
) {
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

    float partial_energies[NUM_OF_THREADS_PER_BLOCK]; // FIXME: remove it?

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
        best_genotype[i] = in_out_genotype[i]; // TODO: add 2dimension
    }

    // Initializing best_energy
    best_energy = INFINITY; // TODO: check correctness

#ifdef ADADELTA_AUTOSTOP
	float rho = 1.0f;
	uint cons_succ = 0;
	uint cons_fail = 0;
#endif


    // TODO: add other blocks for ADADELTA autostop

    // Perfoming ADADELTA iterations

    // The termination criteria is based on a maximum number of iterations,
    // and the minimum step size allowed for single-precision FP numbers
    // (IEEE-754 single float has a precision of about 6 decimal digits)
    do {






#ifdef ADADELTA_AUTOSTOP
	} while ((iteration_cnt < DockConst_max_num_of_iters) && (rho > 0.01f));
#else
    } while (iteration_cnt < DockConst_max_num_of_iters);
#endif

}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
