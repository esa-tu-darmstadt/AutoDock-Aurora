#include "auxiliary.h"

void energy_and_gradient (
	const	uint 				DockConst_pop_size,

			float*				local_coords_x,
			float*				local_coords_y,
			float* 				local_coords_z,

/*
			float				local_coords_x[][MAX_POPSIZE],
			float				local_coords_y[][MAX_POPSIZE],
			float 				local_coords_z[][MAX_POPSIZE],
*/
			float*				gradient_inter_x,
			float*				gradient_inter_y,
			float*				gradient_inter_z,
			float*				gradient_intra_x,
			float*				gradient_intra_y,
			float*				gradient_intra_z,
	// ie
	const 	float* 	restrict	IE_Fgrids,
	const 	float*	restrict	IA_IE_atom_charges,
	const	int*	restrict	IA_IE_atom_types,
			uchar 				DockConst_xsz,
			uchar 				DockConst_ysz,
			uchar 				DockConst_zsz,
			uchar 				DockConst_num_of_atoms,	// TODO: FIXME
			float 				DockConst_gridsize_x_minus1,
			float 				DockConst_gridsize_y_minus1,
			float 				DockConst_gridsize_z_minus1,
			uint 				Host_mul_tmp2,
			uint 				Host_mul_tmp3,
	// ia
	const 	int*	restrict	IA_intraE_contributors,
	const	float*	restrict	IA_reqm,
	const 	float*	restrict	IA_reqm_hbond,
	const 	uint*	restrict	IA_atom1_types_reqm,
	const	uint*	restrict	IA_atom2_types_reqm,
	const	float*	restrict	IA_VWpars_AC,
	const	float*	restrict	IA_VWpars_BD,
	const	float*	restrict	IA_dspars_S,
	const 	float*	restrict	IA_dspars_V,
			float				DockConst_smooth,
			uint 				DockConst_num_of_intraE_contributors,
			float 				DockConst_grid_spacing,
			uint 				DockConst_num_of_atypes,
			float				DockConst_coeff_elec,
			float				DockConst_qasp,
			float				DockConst_coeff_desolv,
    // gradients
    const   int*    restrict    GRAD_rotbonds,
    const   int*    restrict    GRAD_rotbonds_atoms,
    const   int*    restrict    GRAD_num_rotating_atoms_per_rotbond,

    const   float*  restrict    GRAD_angle,
    const   float*  restrict    GRAD_dependence_on_theta,
    const   float*  restrict    GRAD_dependence_on_rotangle
) {

	for (uint atom_id = 0; atom_id < DockConst_num_of_atoms; atom_id++) {

		float x = local_coords_x[atom_id];
		float y = local_coords_y[atom_id];
		float z = local_coords_z[atom_id];

		float partialE1;
		float partialE2;
		float partialE3;

		// If atom is outside of the grid
		if ((x < 0.0f) || (x >= DockConst_gridsize_x_minus1) ||
			(y < 0.0f) || (y >= DockConst_gridsize_y_minus1) ||
			(z < 0.0f) || (z >= DockConst_gridsize_z_minus1)) {

			x = x - 0.5f * DockConst_gridsize_x_minus1; // TODO: in ADGPU is just DockConst_gridsize_x
			y = y - 0.5f * DockConst_gridsize_y_minus1;
			z = z - 0.5f * DockConst_gridsize_z_minus1;

			partialE1 = 16777216.0f;
			partialE2 = 0.0f;
			partialE3 + 0.0f;

			gradient_inter_x[atom_id] = 16777216.0f;
			gradient_inter_y[atom_id] = 16777216.0f;
			gradient_inter_z[atom_id] = 16777216.0f;
		}
		else {
			float x_low = floorf(x);
			float y_low = floorf(y);
			float z_low = floorf(z);

			int ix = (int)x_low;
			int iy = (int)y_low;
			int iz = (int)z_low;

			float dx = x - x_low;
			float dy = y - y_low;
			float dz = z - z_low;

			float omdx = 1.0f - dx;
			float omdy = 1.0f - dy;
			float omdz = 1.0f - dz;

			// Calculating interpolation weights
			float weight000, weight001, weight010;
			float weight011, weight100, weight101;
			float weight110, weight111;
			weight000 = omdx * omdy * omdz;
			weight100 = dx * omdy * omdz;
			weight010 = omdx * dy * omdz;
			weight110 = dx * dy * omdz;
			weight001 = omdx * omdy * dz;
			weight101 = dx * omdy * dz;
			weight011 = omdx * dy * dz;
			weight111 = dx * dy * dz;



		}

	}

}
