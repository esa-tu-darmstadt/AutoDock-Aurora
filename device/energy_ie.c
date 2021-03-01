#include "auxiliary.h"
#include <string.h>

// --------------------------------------------------------------------------
// Calculates the intermolecular energy of a ligand given by
// ligand xyz-positions, and a receptor represented as a grid. 
// The grid point values must be stored at the location which starts at GlobFgrids. 
// If an atom is remains outside the grid, 
// a very high value will be added to the current energy as a penalty. 
// --------------------------------------------------------------------------
void energy_ie (
	const	float*	restrict	IE_Fgrids,
	const	float*	restrict	IA_IE_atom_charges,
	const	int*	restrict	IA_IE_atom_types,
			uchar				DockConst_g1,
			uint				DockConst_g2,
			uint				DockConst_g3,
			uchar				DockConst_num_of_atoms,
			float				DockConst_gridsize_x_minus1,
			float				DockConst_gridsize_y_minus1,
			float				DockConst_gridsize_z_minus1,
			uint				Host_mul_tmp2,
			uint				Host_mul_tmp3,

	        const	uint				DockConst_pop_size,
			float*				final_interE,

			float		local_coords_x[][MAX_POPSIZE],
			float		local_coords_y[][MAX_POPSIZE],
			float		local_coords_z[][MAX_POPSIZE]
)
{
#if defined (PRINT_ALL_IE)
	printf("\n");
	printf("Starting <inter-molecular calculation> ... \n");
	printf("\n");
	printf("\t%-40s %u\n", "DockConst_g1: ", 				DockConst_g1);
	printf("\t%-40s %u\n", "DockConst_g2: ",				DockConst_g2);
	printf("\t%-40s %u\n", "DockConst_g3: ",        		DockConst_g3);
	printf("\t%-40s %u\n", "DockConst_num_of_atoms: ",      DockConst_num_of_atoms);
	printf("\t%-40s %f\n", "DockConst_gridsize_x_minus1: ",	DockConst_gridsize_x_minus1);
	printf("\t%-40s %f\n", "DockConst_gridsize_y_minus1: ", DockConst_gridsize_y_minus1);
	printf("\t%-40s %f\n", "DockConst_gridsize_z_minus1: ", DockConst_gridsize_z_minus1);
	printf("\t%-40s %u\n", "Host_mul_tmp2: ",				Host_mul_tmp2);
	printf("\t%-40s %u\n", "Host_mul_tmp3: ",				Host_mul_tmp3);
#endif

	const float* IE_Fgrids_2 = &IE_Fgrids[Host_mul_tmp2];
	const float* IE_Fgrids_3 = &IE_Fgrids[Host_mul_tmp3];

#if defined (ENABLE_TRACE)
	ftrace_region_begin("IE_MAIN_LOOP");
#endif

#ifdef __clang__
        memset((void *)final_interE, 0, DockConst_pop_size * sizeof(float));
#else
	for (int j = 0; j < DockConst_pop_size; j++) {
		final_interE[j] = 0.0f;
	}
#endif

	// For each ligand atom
#pragma _NEC novector
	for (int atom1_id = 0; atom1_id < DockConst_num_of_atoms; atom1_id++) {
		int atom1_typeid = IA_IE_atom_types[atom1_id];
		float q = IA_IE_atom_charges[atom1_id];

#pragma _NEC vovertake
#pragma _NEC advance_gather
#pragma _NEC gather_reorder
		for (int j = 0; j < DockConst_pop_size; j++) {
     
			float x = local_coords_x[atom1_id][j];
			float y = local_coords_y[atom1_id][j];
			float z = local_coords_z[atom1_id][j];

			float partialE1;
			float partialE2;
			float partialE3;

			// If the atom is outside of the grid
			if ((x < 0.0f) || (x >= DockConst_gridsize_x_minus1) || 
			    (y < 0.0f) || (y >= DockConst_gridsize_y_minus1) ||
			    (z < 0.0f) || (z >= DockConst_gridsize_z_minus1))	{

				// Penalty is 2^24 for each atom outside the grid
				/*
				  interE += 16777216.0f; 
				*/
				partialE1 = 16777216.0f; 
				partialE2 = 0.0f;
				partialE3 = 0.0f;
			}
			else
			{
				int x_low  = (int) floorf(x);
				int y_low  = (int) floorf(y);
				int z_low  = (int) floorf(z);
				int x_high = (int) esa_ceil(x);
				int y_high = (int) esa_ceil(y);
				int z_high = (int) esa_ceil(z);

				float dx = x - x_low;
				float dy = y - y_low;
				float dz = z - z_low;

				// Calculates the weights for trilinear interpolation
				// based on the location of the point inside
				//float weights [2][2][2];
                                float weight000, weight001, weight010;
                                float weight011, weight100, weight101;
                                float weight110, weight111;
				weight000 = (1.0f-dx)*(1.0f-dy)*(1.0f-dz);
				weight100 = dx*(1.0f-dy)*(1.0f-dz);
				weight010 = (1.0f-dx)*dy*(1.0f-dz);
				weight110 = dx*dy*(1.0f-dz);
				weight001 = (1.0f-dx)*(1.0f-dy)*dz;
				weight101 = dx*(1.0f-dy)*dz;
				weight011 = (1.0f-dx)*dy*dz;
				weight111 = dx*dy*dz;

#if defined (PRINT_ALL)
				printf("\n\nPartial results for atom with id %i:\n", atom1_id);
				printf("x_low = %d, x_high = %d, x_frac = %f\n", x_low, x_high, dx);
				printf("y_low = %d, y_high = %d, y_frac = %f\n", y_low, y_high, dy);
				printf("z_low = %d, z_high = %d, z_frac = %f\n\n", z_low, z_high, dz);
				printf("coeff(0,0,0) = %f\n", weight000);
				printf("coeff(1,0,0) = %f\n", weight100);
				printf("coeff(0,1,0) = %f\n", weight010);
				printf("coeff(1,1,0) = %f\n", weight110);
				printf("coeff(0,0,1) = %f\n", weight001);
				printf("coeff(1,0,1) = %f\n", weight101);
				printf("coeff(0,1,1) = %f\n", weight011);
				printf("coeff(1,1,1) = %f\n", weight111);
#endif

				// Added temporal variables
				uint cube_000, cube_100, cube_010, cube_110, cube_001, cube_101, cube_011, cube_111;

				uint ylow_times_g1  = y_low  * DockConst_g1;
				uint yhigh_times_g1 = y_high * DockConst_g1;
				uint zlow_times_g2  = z_low  * DockConst_g2;
				uint zhigh_times_g2 = z_high * DockConst_g2;

				cube_000 = x_low  + ylow_times_g1  + zlow_times_g2;
				cube_100 = x_high + ylow_times_g1  + zlow_times_g2;
				cube_010 = x_low  + yhigh_times_g1 + zlow_times_g2;
				cube_110 = x_high + yhigh_times_g1 + zlow_times_g2;
				cube_001 = x_low  + ylow_times_g1  + zhigh_times_g2;
				cube_101 = x_high + ylow_times_g1  + zhigh_times_g2;
				cube_011 = x_low  + yhigh_times_g1 + zhigh_times_g2;
				cube_111 = x_high + yhigh_times_g1 + zhigh_times_g2;

				uint mul_tmp = atom1_typeid * DockConst_g3;

				// Energy contribution of the current grid type
				//float cube [2][2][2];
                                float cub000, cub001, cub010;
                                float cub011, cub100, cub101;
                                float cub110, cub111;
				cub000 = IE_Fgrids[cube_000 + mul_tmp];
				cub100 = IE_Fgrids[cube_100 + mul_tmp];
				cub010 = IE_Fgrids[cube_010 + mul_tmp];
				cub110 = IE_Fgrids[cube_110 + mul_tmp];
				cub001 = IE_Fgrids[cube_001 + mul_tmp];
				cub101 = IE_Fgrids[cube_101 + mul_tmp];
				cub011 = IE_Fgrids[cube_011 + mul_tmp];
				cub111 = IE_Fgrids[cube_111 + mul_tmp];

#if defined (PRINT_ALL)
				printf("Interpolation of van der Waals map:\n");
				printf("cube(0,0,0) = %f\n", cub000);
				printf("cube(1,0,0) = %f\n", cub100);
				printf("cube(0,1,0) = %f\n", cub010);
				printf("cube(1,1,0) = %f\n", cub110);
				printf("cube(0,0,1) = %f\n", cub001);
				printf("cube(1,0,1) = %f\n", cub101);
				printf("cube(0,1,1) = %f\n", cub011);
				printf("cube(1,1,1) = %f\n", cub111);
#endif

				//partialE1 = TRILININTERPOL(cube, weights);
                                partialE1 =
					cub000 * weight000 +
					cub100 * weight100 +
					cub010 * weight010 +
					cub110 * weight110 +
					cub001 * weight001 +
					cub101 * weight101 +
					cub011 * weight011 +
					cub111 * weight111;

#if defined (PRINT_ALL)
				printf("interpolated value = %f\n\n", TRILININTERPOL(cube, weights));
#endif

				// Energy contribution of the electrostatic grid
				cub000 = IE_Fgrids_2[cube_000];
				cub100 = IE_Fgrids_2[cube_100];
				cub010 = IE_Fgrids_2[cube_010];
				cub110 = IE_Fgrids_2[cube_110];
				cub001 = IE_Fgrids_2[cube_001];
				cub101 = IE_Fgrids_2[cube_101];
				cub011 = IE_Fgrids_2[cube_011];
				cub111 = IE_Fgrids_2[cube_111];

#if defined (PRINT_ALL)
				printf("Interpolation of electrostatic map:\n");
				printf("cube(0,0,0) = %f\n", cub000);
				printf("cube(1,0,0) = %f\n", cub100);
				printf("cube(0,1,0) = %f\n", cub010);
				printf("cube(1,1,0) = %f\n", cub110);
				printf("cube(0,0,1) = %f\n", cub001);
				printf("cube(1,0,1) = %f\n", cub101);
				printf("cube(0,1,1) = %f\n", cub011);
				printf("cube(1,1,1) = %f\n", cub111);
#endif

				//partialE2 = q * TRILININTERPOL(cube, weights);
                                partialE2 = q * (
					cub000 * weight000 +
					cub100 * weight100 +
					cub010 * weight010 +
					cub110 * weight110 +
					cub001 * weight001 +
					cub101 * weight101 +
					cub011 * weight011 +
					cub111 * weight111);

#if defined (PRINT_ALL)
				printf("interpolated value = %f, multiplied by q = %f\n\n", TRILININTERPOL(cube, weights), q*TRILININTERPOL(cube, weights));
#endif

				// Energy contribution of the desolvation grid
				cub000 = IE_Fgrids_3[cube_000];
				cub100 = IE_Fgrids_3[cube_100];
				cub010 = IE_Fgrids_3[cube_010];
				cub110 = IE_Fgrids_3[cube_110];
				cub001 = IE_Fgrids_3[cube_001];
				cub101 = IE_Fgrids_3[cube_101];
				cub011 = IE_Fgrids_3[cube_011];
				cub111 = IE_Fgrids_3[cube_111];

#if defined (PRINT_ALL)
				printf("Interpolation of desolvation map:\n");
				printf("cube(0,0,0) = %f\n", cub000);
				printf("cube(1,0,0) = %f\n", cub100);
				printf("cube(0,1,0) = %f\n", cub010);
				printf("cube(1,1,0) = %f\n", cub110);
				printf("cube(0,0,1) = %f\n", cub001);
				printf("cube(1,0,1) = %f\n", cub101);
				printf("cube(0,1,1) = %f\n", cub011);
				printf("cube(1,1,1) = %f\n", cub111);
#endif

				//partialE3 = fabsf(q) * TRILININTERPOL(cube, weights);
                                partialE3 = fabsf(q) * (
					cub000 * weight000 +
					cub100 * weight100 +
					cub010 * weight010 +
					cub110 * weight110 +
					cub001 * weight001 +
					cub101 * weight101 +
					cub011 * weight011 +
					cub111 * weight111);

#if defined (PRINT_ALL)
				printf("interpolated value = %f, multiplied by abs(q) = %f\n\n", TRILININTERPOL(cube, weights), esa_fabs(q) * trilin_interpol(cube, weights));
				printf("Current value of intermolecular energy = %f\n\n\n", interE);
#endif
			}

			final_interE[j] += partialE1 + partialE2 + partialE3;
	
		} // End j Loop (over individuals)

	} // End of atom1_id for-loop

#if defined (ENABLE_TRACE)
	ftrace_region_end("IE_MAIN_LOOP");
#endif

	// --------------------------------------------------------------

#if defined (PRINT_ALL_IE) 
	printf("\n");
	printf("Finishing <inter-molecular calculation>\n");
	printf("\n");
#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

