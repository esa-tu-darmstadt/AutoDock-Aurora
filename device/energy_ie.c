#include "auxiliary.c"
#include "math.h"

// --------------------------------------------------------------------------
// InterE calculates the intermolecular energy of a ligand given by 
// ligand xyz-positions, and a receptor represented as a grid. 
// The grid point values must be stored at the location which starts at GlobFgrids. 
// If an atom is remains outside the grid, 
// a very high value will be added to the current energy as a penalty. 
// Originally from: processligand.c
// --------------------------------------------------------------------------
void energy_ie (
	const	float*	restrict	GlobFgrids,
	const	float*	restrict	KerConstStatic_atom_charges_const,
	const	char*	restrict	KerConstStatic_atom_types_const,
			uchar				DockConst_g1,
			uint				DockConst_g2,
			uint				DockConst_g3,
			uchar				DockConst_num_of_atoms,
			float				DockConst_gridsize_x_minus1,
			float				DockConst_gridsize_y_minus1,
			float				DockConst_gridsize_z_minus1,
			uint				Host_mul_tmp2,
			uint				Host_mul_tmp3,
			float*				final_interE,

			float*	restrict	local_coords_x,
			float*	restrict	local_coords_y,
			float*	restrict	local_coords_z
)
{
	const float* GlobFgrids2 = & GlobFgrids [Host_mul_tmp2];
	const float* GlobFgrids3 = & GlobFgrids [Host_mul_tmp3];

	#if defined (DEBUG_ACTIVE_KERNEL)
	if (active == 0x00) {printf("	%-20s: %s\n", "Krnl_InterE", "must be disabled");}
	#endif

	float interE = 0.0f;

	// For each ligand atom
	for (uchar atom1_id=0; atom1_id<DockConst_num_of_atoms; atom1_id++)
	{
		char atom1_typeid = KerConstStatic_atom_types_const [atom1_id];

		float x = local_coords_x[atom1_id];
		float y = local_coords_y[atom1_id];
		float z = local_coords_z[atom1_id];
		float q = KerConstStatic_atom_charges_const [atom1_id];

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
			int x_low  = (int) floor(x);
			int y_low  = (int) floor(y);
			int z_low  = (int) floor(z);
			int x_high = (int) ceil(x);	 
			int y_high = (int) ceil(y);
			int z_high = (int) ceil(z);

			float dx = x - x_low; 
			float dy = y - y_low; 
			float dz = z - z_low;

			// Calculates the weights for trilinear interpolation
			// based on the location of the point inside
			float weights [2][2][2];
			weights [0][0][0] = (1-dx)*(1-dy)*(1-dz);
			weights [1][0][0] = dx*(1-dy)*(1-dz);
			weights [0][1][0] = (1-dx)*dy*(1-dz);
			weights [1][1][0] = dx*dy*(1-dz);
			weights [0][0][1] = (1-dx)*(1-dy)*dz;
			weights [1][0][1] = dx*(1-dy)*dz;
			weights [0][1][1] = (1-dx)*dy*dz;
			weights [1][1][1] = dx*dy*dz;

			#if defined (DEBUG_KRNL_INTERE)
			printf("\n\nPartial results for atom with id %i:\n", atom1_id);
			printf("x_low = %d, x_high = %d, x_frac = %f\n", x_low, x_high, dx);
			printf("y_low = %d, y_high = %d, y_frac = %f\n", y_low, y_high, dy);
			printf("z_low = %d, z_high = %d, z_frac = %f\n\n", z_low, z_high, dz);
			printf("coeff(0,0,0) = %f\n", weights [0][0][0]);
			printf("coeff(1,0,0) = %f\n", weights [1][0][0]);
			printf("coeff(0,1,0) = %f\n", weights [0][1][0]);
			printf("coeff(1,1,0) = %f\n", weights [1][1][0]);
			printf("coeff(0,0,1) = %f\n", weights [0][0][1]);
			printf("coeff(1,0,1) = %f\n", weights [1][0][1]);
			printf("coeff(0,1,1) = %f\n", weights [0][1][1]);
			printf("coeff(1,1,1) = %f\n", weights [1][1][1]);
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
			float cube [2][2][2];
	        cube [0][0][0] = GlobFgrids[cube_000 + mul_tmp];
        	cube [1][0][0] = GlobFgrids[cube_100 + mul_tmp];
        	cube [0][1][0] = GlobFgrids[cube_010 + mul_tmp];
        	cube [1][1][0] = GlobFgrids[cube_110 + mul_tmp];
        	cube [0][0][1] = GlobFgrids[cube_001 + mul_tmp];
        	cube [1][0][1] = GlobFgrids[cube_101 + mul_tmp];
        	cube [0][1][1] = GlobFgrids[cube_011 + mul_tmp];
        	cube [1][1][1] = GlobFgrids[cube_111 + mul_tmp];
		
			#if defined (DEBUG_KRNL_INTERE)
			printf("Interpolation of van der Waals map:\n");
			printf("cube(0,0,0) = %f\n", cube [0][0][0]);
			printf("cube(1,0,0) = %f\n", cube [1][0][0]);
			printf("cube(0,1,0) = %f\n", cube [0][1][0]);
			printf("cube(1,1,0) = %f\n", cube [1][1][0]);
			printf("cube(0,0,1) = %f\n", cube [0][0][1]);
			printf("cube(1,0,1) = %f\n", cube [1][0][1]);
			printf("cube(0,1,1) = %f\n", cube [0][1][1]);
			printf("cube(1,1,1) = %f\n", cube [1][1][1]);
			#endif

			/*partialE1 = TRILININTERPOL(cube, weights);*/
			partialE1 = cube[0][0][0] * weights[0][0][0] +
				    cube[1][0][0] * weights[1][0][0] +
				    cube[0][1][0] * weights[0][1][0] +
				    cube[1][1][0] * weights[1][1][0] + 
				    cube[0][0][1] * weights[0][0][1] +
				    cube[1][0][1] * weights[1][0][1] + 
				    cube[0][1][1] * weights[0][1][1] +
				    cube[1][1][1] * weights[1][1][1];

			#if defined (DEBUG_KRNL_INTERE)
			printf("interpolated value = %f\n\n", TRILININTERPOL(cube, weights));
			#endif

			// Energy contribution of the electrostatic grid
			cube [0][0][0] = GlobFgrids2[cube_000] /*GlobFgrids [Host_mul_tmp2 + cube_000]*/;
                        cube [1][0][0] = GlobFgrids2[cube_100] /*GlobFgrids [Host_mul_tmp2 + cube_100]*/;
                        cube [0][1][0] = GlobFgrids2[cube_010] /*GlobFgrids [Host_mul_tmp2 + cube_010]*/;
                        cube [1][1][0] = GlobFgrids2[cube_110] /*GlobFgrids [Host_mul_tmp2 + cube_110]*/;
                        cube [0][0][1] = GlobFgrids2[cube_001] /*GlobFgrids [Host_mul_tmp2 + cube_001]*/;
                        cube [1][0][1] = GlobFgrids2[cube_101] /*GlobFgrids [Host_mul_tmp2 + cube_101]*/;
                        cube [0][1][1] = GlobFgrids2[cube_011] /*GlobFgrids [Host_mul_tmp2 + cube_011]*/;
                        cube [1][1][1] = GlobFgrids2[cube_111] /*GlobFgrids [Host_mul_tmp2 + cube_111]*/;

			#if defined (DEBUG_KRNL_INTERE)
			printf("Interpolation of electrostatic map:\n");
			printf("cube(0,0,0) = %f\n", cube [0][0][0]);
			printf("cube(1,0,0) = %f\n", cube [1][0][0]);
			printf("cube(0,1,0) = %f\n", cube [0][1][0]);
			printf("cube(1,1,0) = %f\n", cube [1][1][0]);
			printf("cube(0,0,1) = %f\n", cube [0][0][1]);
			printf("cube(1,0,1) = %f\n", cube [1][0][1]);
			printf("cube(0,1,1) = %f\n", cube [0][1][1]);
			printf("cube(1,1,1) = %f\n", cube [1][1][1]);
			#endif

			/*partialE2 = q * TRILININTERPOL(cube, weights);*/
			partialE2 = q * (cube[0][0][0] * weights[0][0][0] +
				    	 cube[1][0][0] * weights[1][0][0] +
				    	 cube[0][1][0] * weights[0][1][0] +
				    	 cube[1][1][0] * weights[1][1][0] + 
				    	 cube[0][0][1] * weights[0][0][1] +
				    	 cube[1][0][1] * weights[1][0][1] + 
				    	 cube[0][1][1] * weights[0][1][1] +
				    	 cube[1][1][1] * weights[1][1][1]);
		
			#if defined (DEBUG_KRNL_INTERE)
			printf("interpolated value = %f, multiplied by q = %f\n\n", TRILININTERPOL(cube, weights), q*TRILININTERPOL(cube, weights));
			#endif

			// Energy contribution of the desolvation grid
			cube [0][0][0] = GlobFgrids3[cube_000] /*GlobFgrids [Host_mul_tmp3 + cube_000]*/;
                        cube [1][0][0] = GlobFgrids3[cube_100] /*GlobFgrids [Host_mul_tmp3 + cube_100]*/;
                        cube [0][1][0] = GlobFgrids3[cube_010] /*GlobFgrids [Host_mul_tmp3 + cube_010]*/;
                        cube [1][1][0] = GlobFgrids3[cube_110] /*GlobFgrids [Host_mul_tmp3 + cube_110]*/;
                        cube [0][0][1] = GlobFgrids3[cube_001] /*GlobFgrids [Host_mul_tmp3 + cube_001]*/;
                        cube [1][0][1] = GlobFgrids3[cube_101] /*GlobFgrids [Host_mul_tmp3 + cube_101]*/;
                        cube [0][1][1] = GlobFgrids3[cube_011] /*GlobFgrids [Host_mul_tmp3 + cube_011]*/;
                        cube [1][1][1] = GlobFgrids3[cube_111] /*GlobFgrids [Host_mul_tmp3 + cube_111]*/;

			#if defined (DEBUG_KRNL_INTERE)
			printf("Interpolation of desolvation map:\n");
			printf("cube(0,0,0) = %f\n", cube [0][0][0]);
			printf("cube(1,0,0) = %f\n", cube [1][0][0]);
			printf("cube(0,1,0) = %f\n", cube [0][1][0]);
			printf("cube(1,1,0) = %f\n", cube [1][1][0]);
			printf("cube(0,0,1) = %f\n", cube [0][0][1]);
			printf("cube(1,0,1) = %f\n", cube [1][0][1]);
			printf("cube(0,1,1) = %f\n", cube [0][1][1]);
			printf("cube(1,1,1) = %f\n", cube [1][1][1]);
			#endif

			/*partialE3 = fabs(q) * TRILININTERPOL(cube, weights);*/
			partialE3 = fabs(q) * (cube[0][0][0] * weights[0][0][0] +
				    	       cube[1][0][0] * weights[1][0][0] +
				    	       cube[0][1][0] * weights[0][1][0] +
				    	       cube[1][1][0] * weights[1][1][0] + 
				    	       cube[0][0][1] * weights[0][0][1] +
				    	       cube[1][0][1] * weights[1][0][1] + 
				    	       cube[0][1][1] * weights[0][1][1] +
				    	       cube[1][1][1] * weights[1][1][1]);

			#if defined (DEBUG_KRNL_INTERE)
			printf("interpolated value = %f, multiplied by abs(q) = %f\n\n", TRILININTERPOL(cube, weights), fabs(q) * trilin_interpol(cube, weights));
			printf("Current value of intermolecular energy = %f\n\n\n", interE);
			#endif
		}

		interE += partialE1 + partialE2 + partialE3;
	} // End of atom1_id for-loop

	// --------------------------------------------------------------
	// Send intermolecular energy to chanel
	// --------------------------------------------------------------
	*final_interE = interE;

	// --------------------------------------------------------------

	#if defined (DEBUG_KRNL_INTERE)
	printf("AFTER Out INTERE CHANNEL\n");
	#endif
 	
	#if defined (DEBUG_ACTIVE_KERNEL)
	printf("	%-20s: %s\n", "Krnl_InterE", "disabled");
	#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

