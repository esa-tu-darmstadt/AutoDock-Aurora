#include "../defines.h"

// Constant struct
typedef struct
{
       float atom_charges_const[MAX_NUM_OF_ATOMS];
       char  atom_types_const  [MAX_NUM_OF_ATOMS];
       char  intraE_contributors_const[3*MAX_INTRAE_CONTRIBUTORS];
       float VWpars_AC_const   [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
       float VWpars_BD_const   [MAX_NUM_OF_ATYPES*MAX_NUM_OF_ATYPES];
       float dspars_S_const    [MAX_NUM_OF_ATYPES];
       float dspars_V_const    [MAX_NUM_OF_ATYPES];
       int   rotlist_const     [MAX_NUM_OF_ROTATIONS];
       float ref_coords_x_const[MAX_NUM_OF_ATOMS];
       float ref_coords_y_const[MAX_NUM_OF_ATOMS];
       float ref_coords_z_const[MAX_NUM_OF_ATOMS];
       float rotbonds_moving_vectors_const[3*MAX_NUM_OF_ROTBONDS];
       float rotbonds_unit_vectors_const  [3*MAX_NUM_OF_ROTBONDS];
       float ref_orientation_quats_const  [4*MAX_NUM_OF_RUNS];
} kernelconstant;

// All related pragmas are in defines.h (accesible by host and device code)

// -------------------------------------------------------
//
// -------------------------------------------------------
void calc_energy(int    dockpars_rotbondlist_length,
		 char   dockpars_num_of_atoms,
		 char   dockpars_gridsize_x,
		 char   dockpars_gridsize_y,
                 char   dockpars_gridsize_z,
  __global const float* restrict dockpars_fgrids, // cannot be allocated in __constant (too large)
	         char   dockpars_num_of_atypes,
		 int    dockpars_num_of_intraE_contributors,
		 float  dockpars_grid_spacing,
		 float  dockpars_coeff_elec,
                 float  dockpars_qasp,
		 float  dockpars_coeff_desolv,
         __local float* genotype, 
	 __local float* energy,
	 __local int*   run_id,
   	 // Altera doesn't allow local var outside kernels
 	 // so this local vars are passed from a kernel
	 __local float* calc_coords_x,
	 __local float* calc_coords_y,
	 __local float* calc_coords_z,
	 __local float* partial_energies,

__global kernelconstant* restrict KerConst

/*
      __constant float* restrict atom_charges_const,
      __constant char*  restrict atom_types_const,
      __constant char*  restrict intraE_contributors_const,
      __constant float* restrict VWpars_AC_const,
      __constant float* restrict VWpars_BD_const,
      __constant float* restrict dspars_S_const,
      __constant float* restrict dspars_V_const,
      __constant int*   restrict rotlist_const,
      __constant float* restrict ref_coords_x_const,
      __constant float* restrict ref_coords_y_const,
      __constant float* restrict ref_coords_z_const,
      __constant float* restrict rotbonds_moving_vectors_const,
      __constant float* restrict rotbonds_unit_vectors_const,
      __constant float* restrict ref_orientation_quats_const
*/
)

//The GPU device function calculates the energy of the entity described by genotype, dockpars and the liganddata
//arrays in constant memory and returns it in the energy parameter. The parameter run_id has to be equal to the ID
//of the run whose population includes the current entity (which can be determined with blockIdx.x), since this
//determines which reference orientation should be used.
{

	int contributor_counter;
	char atom1_id, atom2_id, atom1_typeid, atom2_typeid;

	// Name changed to distance_leo to avoid
	// errors as "distance" is the name of OpenCL function
	//float subx, suby, subz, distance;
	float subx, suby, subz, distance_leo;

	float x, y, z, dx, dy, dz, q;
	float cube[2][2][2];
	float weights[2][2][2];
	int x_low, x_high, y_low, y_high, z_low, z_high;

	float phi, theta, genrotangle, rotation_angle, sin_angle;
	float genrot_unitvec[3], rotation_unitvec[3], rotation_movingvec[3];
	int rotation_counter, rotation_list_element;
	float atom_to_rotate[3];
	int atom_id, rotbond_id;
	float quatrot_left_x, quatrot_left_y, quatrot_left_z, quatrot_left_q;
	float quatrot_temp_x, quatrot_temp_y, quatrot_temp_z, quatrot_temp_q;

	// Altera doesn't allow local var outside kernels
	// so this local vars are passed from a kernel
	//__local float calc_coords_x[MAX_NUM_OF_ATOMS];
	//__local float calc_coords_y[MAX_NUM_OF_ATOMS];
	//__local float calc_coords_z[MAX_NUM_OF_ATOMS];
	//__local float partial_energies[NUM_OF_THREADS_PER_BLOCK];

	partial_energies[get_local_id(0)] = 0.0f;

	//CALCULATE CONFORMATION

	//calculate vectors for general rotation
	phi         = genotype[3]*DEG_TO_RAD;
	theta       = genotype[4]*DEG_TO_RAD;
	genrotangle = genotype[5]*DEG_TO_RAD;

#if defined (IMPROVE_GRID)
	#if defined (NATIVE_PRECISION)
	sin_angle = native_sin(theta);
	genrot_unitvec [0] = sin_angle*native_cos(phi);
	genrot_unitvec [1] = sin_angle*native_sin(phi);
	genrot_unitvec [2] = native_cos(theta);
	#elif defined (HALF_PRECISION)
	sin_angle = half_sin(theta);
	genrot_unitvec [0] = sin_angle*half_cos(phi);
	genrot_unitvec [1] = sin_angle*half_sin(phi);
	genrot_unitvec [2] = half_cos(theta);
	#else	// Full precision
	sin_angle = sin(theta);
	genrot_unitvec [0] = sin_angle*cos(phi);
	genrot_unitvec [1] = sin_angle*sin(phi);
	genrot_unitvec [2] = cos(theta);
	#endif

	// INTERMOLECULAR for-loop (intermediate results)
	// It stores a product of two chars
	unsigned int mul_tmp;

	unsigned char g1 = dockpars_gridsize_x;
	unsigned int  g2 = dockpars_gridsize_x * dockpars_gridsize_y;
        unsigned int  g3 = dockpars_gridsize_x * dockpars_gridsize_y * dockpars_gridsize_z;

	unsigned int ylow_times_g1, yhigh_times_g1;
	unsigned int zlow_times_g2, zhigh_times_g2;

	unsigned int cube_000;
	unsigned int cube_100;
        unsigned int cube_010;
	unsigned int cube_110;
        unsigned int cube_001;
        unsigned int cube_101;
        unsigned int cube_011;
        unsigned int cube_111;
#else
	sin_angle = sin(theta);
	genrot_unitvec [0] = sin_angle*cos(phi);
	genrot_unitvec [1] = sin_angle*sin(phi);
	genrot_unitvec [2] = cos(theta);
#endif


	// ================================================
	// Iterating over elements of rotation list
	// ================================================

	for (rotation_counter = get_local_id(0);
	     rotation_counter < dockpars_rotbondlist_length;
	     rotation_counter+=NUM_OF_THREADS_PER_BLOCK)
	{
		rotation_list_element = KerConst->rotlist_const[rotation_counter];

		if ((rotation_list_element & RLIST_DUMMY_MASK) == 0)	//if not dummy rotation
		{
			atom_id = rotation_list_element & RLIST_ATOMID_MASK;

			//capturing atom coordinates
			if ((rotation_list_element & RLIST_FIRSTROT_MASK) != 0)	//if firts rotation of this atom
			{
				atom_to_rotate[0] = KerConst->ref_coords_x_const[atom_id];
				atom_to_rotate[1] = KerConst->ref_coords_y_const[atom_id];
				atom_to_rotate[2] = KerConst->ref_coords_z_const[atom_id];
			}
			else
			{
				atom_to_rotate[0] = calc_coords_x[atom_id];
				atom_to_rotate[1] = calc_coords_y[atom_id];
				atom_to_rotate[2] = calc_coords_z[atom_id];
			}

			//capturing rotation vectors and angle
			if ((rotation_list_element & RLIST_GENROT_MASK) != 0)	//if general rotation
			{
				rotation_unitvec[0] = genrot_unitvec[0];
				rotation_unitvec[1] = genrot_unitvec[1];
				rotation_unitvec[2] = genrot_unitvec[2];

				rotation_angle = genrotangle;

				rotation_movingvec[0] = genotype[0];
				rotation_movingvec[1] = genotype[1];
				rotation_movingvec[2] = genotype[2];
			}
			else	//if rotating around rotatable bond
			{
				rotbond_id = (rotation_list_element & RLIST_RBONDID_MASK) >> RLIST_RBONDID_SHIFT;
	
				rotation_unitvec[0] = KerConst->rotbonds_unit_vectors_const[3*rotbond_id];
				rotation_unitvec[1] = KerConst->rotbonds_unit_vectors_const[3*rotbond_id+1];
				rotation_unitvec[2] = KerConst->rotbonds_unit_vectors_const[3*rotbond_id+2];
				rotation_angle = genotype[6+rotbond_id]*DEG_TO_RAD;

				rotation_movingvec[0] = KerConst->rotbonds_moving_vectors_const[3*rotbond_id];
				rotation_movingvec[1] = KerConst->rotbonds_moving_vectors_const[3*rotbond_id+1];
				rotation_movingvec[2] = KerConst->rotbonds_moving_vectors_const[3*rotbond_id+2];

				//in addition performing the first movement 
				//which is needed only if rotating around rotatable bond
				atom_to_rotate[0] -= rotation_movingvec[0];
				atom_to_rotate[1] -= rotation_movingvec[1];
				atom_to_rotate[2] -= rotation_movingvec[2];
			}

			//performing rotation

#if defined (NATIVE_PRECISION)		
			rotation_angle = native_divide(rotation_angle,2);
			quatrot_left_q = native_cos(rotation_angle);
			sin_angle = native_sin(rotation_angle);
#elif defined (HALF_PRECISION)
			rotation_angle = half_divide(rotation_angle,2);
			quatrot_left_q = half_cos(rotation_angle);
			sin_angle = half_sin(rotation_angle);
#else	// Full precision
			rotation_angle = rotation_angle/2;
			quatrot_left_q = cos(rotation_angle);
			sin_angle = sin(rotation_angle);
#endif
			quatrot_left_x = sin_angle*rotation_unitvec[0];
			quatrot_left_y = sin_angle*rotation_unitvec[1];
			quatrot_left_z = sin_angle*rotation_unitvec[2];

			if ((rotation_list_element & RLIST_GENROT_MASK) != 0)	//if general rotation, 
										//two rotations should be performed 
										//(multiplying the quaternions)
			{
				//calculating quatrot_left*ref_orientation_quats_const, 
				//which means that reference orientation rotation is the first
				quatrot_temp_q = quatrot_left_q;
				quatrot_temp_x = quatrot_left_x;
				quatrot_temp_y = quatrot_left_y;
				quatrot_temp_z = quatrot_left_z;

				quatrot_left_q = quatrot_temp_q*KerConst->ref_orientation_quats_const[4*(*run_id)]-
						 quatrot_temp_x*KerConst->ref_orientation_quats_const[4*(*run_id)+1]-
						 quatrot_temp_y*KerConst->ref_orientation_quats_const[4*(*run_id)+2]-
						 quatrot_temp_z*KerConst->ref_orientation_quats_const[4*(*run_id)+3];
				quatrot_left_x = quatrot_temp_q*KerConst->ref_orientation_quats_const[4*(*run_id)+1]+
						 KerConst->ref_orientation_quats_const[4*(*run_id)]*quatrot_temp_x+
						 quatrot_temp_y*KerConst->ref_orientation_quats_const[4*(*run_id)+3]-
						 KerConst->ref_orientation_quats_const[4*(*run_id)+2]*quatrot_temp_z;
				quatrot_left_y = quatrot_temp_q*KerConst->ref_orientation_quats_const[4*(*run_id)+2]+
						 KerConst->ref_orientation_quats_const[4*(*run_id)]*quatrot_temp_y+
						 KerConst->ref_orientation_quats_const[4*(*run_id)+1]*quatrot_temp_z-
						 quatrot_temp_x*KerConst->ref_orientation_quats_const[4*(*run_id)+3];
				quatrot_left_z = quatrot_temp_q*KerConst->ref_orientation_quats_const[4*(*run_id)+3]+
						 KerConst->ref_orientation_quats_const[4*(*run_id)]*quatrot_temp_z+
						 quatrot_temp_x*KerConst->ref_orientation_quats_const[4*(*run_id)+2]-
						 KerConst->ref_orientation_quats_const[4*(*run_id)+1]*quatrot_temp_y;

			}

			quatrot_temp_q = 0 -
					 quatrot_left_x*atom_to_rotate [0] -
					 quatrot_left_y*atom_to_rotate [1] -
					 quatrot_left_z*atom_to_rotate [2];
			quatrot_temp_x = quatrot_left_q*atom_to_rotate [0] +
					 quatrot_left_y*atom_to_rotate [2] -
					 quatrot_left_z*atom_to_rotate [1];
			quatrot_temp_y = quatrot_left_q*atom_to_rotate [1] -
					 quatrot_left_x*atom_to_rotate [2] +
					 quatrot_left_z*atom_to_rotate [0];
			quatrot_temp_z = quatrot_left_q*atom_to_rotate [2] +
					 quatrot_left_x*atom_to_rotate [1] -
					 quatrot_left_y*atom_to_rotate [0];

			atom_to_rotate [0] = 0 -
					     quatrot_temp_q*quatrot_left_x +
					     quatrot_temp_x*quatrot_left_q -
					     quatrot_temp_y*quatrot_left_z +
					     quatrot_temp_z*quatrot_left_y;
			atom_to_rotate [1] = 0 -
					     quatrot_temp_q*quatrot_left_y +
					     quatrot_temp_x*quatrot_left_z +
					     quatrot_temp_y*quatrot_left_q -
					     quatrot_temp_z*quatrot_left_x;
			atom_to_rotate [2] = 0 -
					     quatrot_temp_q*quatrot_left_z -
					     quatrot_temp_x*quatrot_left_y +
					     quatrot_temp_y*quatrot_left_x +
					     quatrot_temp_z*quatrot_left_q;

			//performing final movement and storing values
			calc_coords_x[atom_id] = atom_to_rotate [0] + rotation_movingvec[0];
			calc_coords_y[atom_id] = atom_to_rotate [1] + rotation_movingvec[1];
			calc_coords_z[atom_id] = atom_to_rotate [2] + rotation_movingvec[2];

		} // End if-statement not dummy rotation

		barrier(CLK_LOCAL_MEM_FENCE);

	} // End rotation_counter for-loop
	

	// ================================================
	// CALCULATE INTERMOLECULAR ENERGY
	// ================================================

	for (atom1_id = get_local_id(0);
	     atom1_id < dockpars_num_of_atoms;
	     atom1_id+= NUM_OF_THREADS_PER_BLOCK)
	{
		atom1_typeid = KerConst->atom_types_const[atom1_id];
		x = calc_coords_x[atom1_id];
		y = calc_coords_y[atom1_id];
		z = calc_coords_z[atom1_id];
		q = KerConst->atom_charges_const[atom1_id];

		if ((x < 0) || (y < 0) || (z < 0) || (x >= dockpars_gridsize_x-1)
				                  || (y >= dockpars_gridsize_y-1)
						  || (z >= dockpars_gridsize_z-1)){

			partial_energies[get_local_id(0)] += 16777216.0f; //100000.0f;
		}
		else
		{
			//get coordinates
			x_low = (int)floor(x); y_low = (int)floor(y); z_low = (int)floor(z); 
			x_high = (int)ceil(x); y_high = (int)ceil(y); z_high = (int)ceil(z);
			dx = x - x_low; dy = y - y_low; dz = z - z_low;

			//calculate interpolation weights
			weights [0][0][0] = (1-dx)*(1-dy)*(1-dz);
			weights [1][0][0] = dx*(1-dy)*(1-dz);
			weights [0][1][0] = (1-dx)*dy*(1-dz);
			weights [1][1][0] = dx*dy*(1-dz);
			weights [0][0][1] = (1-dx)*(1-dy)*dz;
			weights [1][0][1] = dx*(1-dy)*dz;
			weights [0][1][1] = (1-dx)*dy*dz;
			weights [1][1][1] = dx*dy*dz;

			//capturing affinity values
#if defined (IMPROVE_GRID)
			ylow_times_g1  = y_low*g1;
			yhigh_times_g1 = y_high*g1;
		        zlow_times_g2  = z_low*g2;
			zhigh_times_g2 = z_high*g2;
			
			cube_000 = x_low  + ylow_times_g1  + zlow_times_g2; 		
			cube_100 = x_high + ylow_times_g1  + zlow_times_g2;
			cube_010 = x_low  + yhigh_times_g1 + zlow_times_g2;	
			cube_110 = x_high + yhigh_times_g1 + zlow_times_g2;
			cube_001 = x_low  + ylow_times_g1  + zhigh_times_g2;
			cube_101 = x_high + ylow_times_g1  + zhigh_times_g2;
			cube_011 = x_low  + yhigh_times_g1 + zhigh_times_g2;
			cube_111 = x_high + yhigh_times_g1 + zhigh_times_g2;
			mul_tmp = atom1_typeid*g3;

			cube [0][0][0] = *(dockpars_fgrids + cube_000 + mul_tmp);
			cube [1][0][0] = *(dockpars_fgrids + cube_100 + mul_tmp);
                        cube [0][1][0] = *(dockpars_fgrids + cube_010 + mul_tmp);
                        cube [1][1][0] = *(dockpars_fgrids + cube_110 + mul_tmp);
                        cube [0][0][1] = *(dockpars_fgrids + cube_001 + mul_tmp);
                        cube [1][0][1] = *(dockpars_fgrids + cube_101 + mul_tmp);
                        cube [0][1][1] = *(dockpars_fgrids + cube_011 + mul_tmp);
                        cube [1][1][1] = *(dockpars_fgrids + cube_111 + mul_tmp);

#else
			cube [0][0][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_low, x_low);
			cube [1][0][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_low, x_high);
			cube [0][1][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_high, x_low);
			cube [1][1][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_high, x_high);
			cube [0][0][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_low, x_low);
			cube [1][0][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_low, x_high);
			cube [0][1][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_high, x_low);
			cube [1][1][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_high, x_high);
#endif
			//calculating affinity energy
			partial_energies[get_local_id(0)] += TRILININTERPOL(cube, weights);

			//capturing electrostatic values
			atom1_typeid = dockpars_num_of_atypes;

#if defined (IMPROVE_GRID)
			mul_tmp = atom1_typeid*g3;
			cube [0][0][0] = *(dockpars_fgrids + cube_000 + mul_tmp);
			cube [1][0][0] = *(dockpars_fgrids + cube_100 + mul_tmp);
                        cube [0][1][0] = *(dockpars_fgrids + cube_010 + mul_tmp);
                        cube [1][1][0] = *(dockpars_fgrids + cube_110 + mul_tmp);
                        cube [0][0][1] = *(dockpars_fgrids + cube_001 + mul_tmp);
                        cube [1][0][1] = *(dockpars_fgrids + cube_101 + mul_tmp);
                        cube [0][1][1] = *(dockpars_fgrids + cube_011 + mul_tmp);
                        cube [1][1][1] = *(dockpars_fgrids + cube_111 + mul_tmp);

#else
			cube [0][0][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_low, x_low);
			cube [1][0][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_low, x_high);
			cube [0][1][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_high, x_low);
			cube [1][1][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_high, x_high);
			cube [0][0][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_low, x_low);
			cube [1][0][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_low, x_high);
			cube [0][1][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_high, x_low);
			cube [1][1][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_high, x_high);
#endif
			//calculating electrosatic energy
			partial_energies[get_local_id(0)] += q * TRILININTERPOL(cube, weights);

			//capturing desolvation values
			atom1_typeid = dockpars_num_of_atypes+1;

#if defined (IMPROVE_GRID)
			mul_tmp = atom1_typeid*g3;
			cube [0][0][0] = *(dockpars_fgrids + cube_000 + mul_tmp);
			cube [1][0][0] = *(dockpars_fgrids + cube_100 + mul_tmp);
                        cube [0][1][0] = *(dockpars_fgrids + cube_010 + mul_tmp);
                        cube [1][1][0] = *(dockpars_fgrids + cube_110 + mul_tmp);
                        cube [0][0][1] = *(dockpars_fgrids + cube_001 + mul_tmp);
                        cube [1][0][1] = *(dockpars_fgrids + cube_101 + mul_tmp);
                        cube [0][1][1] = *(dockpars_fgrids + cube_011 + mul_tmp);
                        cube [1][1][1] = *(dockpars_fgrids + cube_111 + mul_tmp);

#else
			cube [0][0][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_low, x_low);
			cube [1][0][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_low, x_high);
			cube [0][1][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_high, x_low);
			cube [1][1][0] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_low, y_high, x_high);
			cube [0][0][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_low, x_low);
			cube [1][0][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_low, x_high);
			cube [0][1][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_high, x_low);
			cube [1][1][1] = GETGRIDVALUE(dockpars_fgrids, dockpars_gridsize_x,
						      dockpars_gridsize_y, dockpars_gridsize_z,
						      atom1_typeid, z_high, y_high, x_high);
#endif
			//calculating desolvation energy
			partial_energies[get_local_id(0)] += fabs(q) * TRILININTERPOL(cube, weights);
		}

	} // End atom1_id for-loop


	// Altera doesn't support power function
	// so this is implemented with multiplications
	// Full precision is used
	float distance_pow_2, distance_pow_4, distance_pow_6, distance_pow_10, distance_pow_12;


	// In paper: intermolecular and internal energy calculation
	// are independent from each other, -> NO BARRIER NEEDED
        // but require different operations, 
	// thus, they can be executed only sequentially on the GPU.

	// ================================================
	// CALCULATE INTRAMOLECULAR ENERGY
	// ================================================

	for (contributor_counter = get_local_id(0);
	     contributor_counter < dockpars_num_of_intraE_contributors;
	     contributor_counter +=NUM_OF_THREADS_PER_BLOCK)
	{

		//getting atom IDs
		atom1_id = KerConst->intraE_contributors_const[3*contributor_counter];
		atom2_id = KerConst->intraE_contributors_const[3*contributor_counter+1];

		//atom1_id = KerConst->intraE_contributors_const[0];
		//atom2_id = KerConst->intraE_contributors_const[1];

		//calculating address of first atom's coordinates
		subx = calc_coords_x[atom1_id];
		suby = calc_coords_y[atom1_id];
		subz = calc_coords_z[atom1_id];

		//calculating address of second atom's coordinates
		subx -= calc_coords_x[atom2_id];
		suby -= calc_coords_y[atom2_id];
		subz -= calc_coords_z[atom2_id];

		//calculating distance (distance_leo)
#if defined (NATIVE_PRECISION)
		distance_leo = native_sqrt(subx*subx + suby*suby + subz*subz)*dockpars_grid_spacing;
#elif defined (HALF_PRECISION)
		distance_leo = half_sqrt(subx*subx + suby*suby + subz*subz)*dockpars_grid_spacing;
#else	// Full precision
		distance_leo = sqrt(subx*subx + suby*suby + subz*subz)*dockpars_grid_spacing;
#endif

		if (distance_leo < 1.0f)
			distance_leo = 1.0f;

		distance_pow_2 = distance_leo*distance_leo;
		distance_pow_4 = distance_pow_2*distance_pow_2;
		distance_pow_6 = distance_pow_2*distance_pow_4;
		distance_pow_10 = distance_pow_4*distance_pow_6;
		distance_pow_12 = distance_pow_6*distance_pow_6;

		//calculating energy contributions
		if ((distance_leo < 8.0f) && (distance_leo < 20.48f))
		{
			//getting type IDs
			atom1_typeid = KerConst->atom_types_const[atom1_id];
			atom2_typeid = KerConst->atom_types_const[atom2_id];

			//calculating van der Waals / hydrogen bond term
#if defined (NATIVE_PRECISION)
			partial_energies[get_local_id(0)] += native_divide(KerConst->VWpars_AC_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid],native_powr(distance_leo,12));
#elif defined (HALF_PRECISION)
			partial_energies[get_local_id(0)] += half_divide(KerConst->VWpars_AC_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid],half_powr(distance_leo,12));
#else	// Full precision // from 79 to 109 sec
			//partial_energies[get_local_id(0)] += KerConst->VWpars_AC_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid]/powr(distance_leo,12);
			partial_energies[get_local_id(0)] += KerConst->VWpars_AC_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid]/distance_pow_12;
#endif

			if (KerConst->intraE_contributors_const[3*contributor_counter+2] == 1)	//H-bond
#if defined (NATIVE_PRECISION)
				partial_energies[get_local_id(0)] -= native_divide(KerConst->VWpars_BD_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid],native_powr(distance_leo,10));
#elif defined (HALF_PRECISION)
				partial_energies[get_local_id(0)] -= half_divide(KerConst->VWpars_BD_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid],half_powr(distance_leo,10));
#else	// Full precision // from 109 to 121 sec
				//partial_energies[get_local_id(0)] -= KerConst->VWpars_BD_const[atom1_typeid*dockpars_num_of_atypes+atom2_typeid]/powr(distance_leo,10);
				partial_energies[get_local_id(0)] -= KerConst->VWpars_BD_const[atom1_typeid*dockpars_num_of_atypes+atom2_typeid]/distance_pow_10;
#endif

			else	//van der Waals
#if defined (NATIVE_PRECISION)
				partial_energies[get_local_id(0)] -= native_divide(KerConst->VWpars_BD_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid],native_powr(distance_leo,6));
#elif defined (HALF_PRECISION)
				partial_energies[get_local_id(0)] -= half_divide(KerConst->VWpars_BD_const[atom1_typeid * dockpars_num_of_atypes+atom2_typeid],half_powr(distance_leo,6));
#else	// Full precision // from 121 to 146 sec
				//partial_energies[get_local_id(0)] -= KerConst->VWpars_BD_const[atom1_typeid*dockpars_num_of_atypes+atom2_typeid]/powr(distance_leo,6);
				partial_energies[get_local_id(0)] -= KerConst->VWpars_BD_const[atom1_typeid*dockpars_num_of_atypes+atom2_typeid]/distance_pow_6;
#endif

			//calculating electrostatic term
#if defined (NATIVE_PRECISION)
                        partial_energies[get_local_id(0)] += native_divide (
                                                             dockpars_coeff_elec * KerConst->atom_charges_const[atom1_id] * KerConst->atom_charges_const[atom2_id],
                                                             distance_leo * (-8.5525f + native_divide(86.9525f,(1.0f + 7.7839f*native_exp(-0.3154f*distance_leo))))
                                                             );
#elif defined (HALF_PRECISION)
                        partial_energies[get_local_id(0)] += half_divide (
                                                             dockpars_coeff_elec * KerConst->atom_charges_const[atom1_id] * KerConst->atom_charges_const[atom2_id],
                                                             distance_leo * (-8.5525f + half_divide(86.9525f,(1.0f + 7.7839f*half_exp(-0.3154f*distance_leo))))
                                                             );
#else	// Full precision // from 146 to 148 sec
			partial_energies[get_local_id(0)] += dockpars_coeff_elec*KerConst->atom_charges_const[atom1_id]*KerConst->atom_charges_const[atom2_id]/
			                                     (distance_leo*(-8.5525f + 86.9525f/(1.0f + 7.7839f*exp(-0.3154f*distance_leo))));
#endif

			//calculating desolvation term
#if defined (NATIVE_PRECISION)
			partial_energies[get_local_id(0)] += ((KerConst->dspars_S_const[atom1_typeid] + 
							       dockpars_qasp*fabs(KerConst->atom_charges_const[atom1_id]))*KerConst->dspars_V_const[atom2_typeid] +
					                      (KerConst->dspars_S_const[atom2_typeid] + 
							       dockpars_qasp*fabs(KerConst->atom_charges_const[atom2_id]))*KerConst->dspars_V_const[atom1_typeid]) *
					                       dockpars_coeff_desolv*native_exp(-distance_leo*native_divide(distance_leo,25.92f));
#elif defined (HALF_PRECISION)
			partial_energies[get_local_id(0)] += ((KerConst->dspars_S_const[atom1_typeid] + 
							       dockpars_qasp*fabs(KerConst->atom_charges_const[atom1_id]))*KerConst->dspars_V_const[atom2_typeid] +
					                      (KerConst->dspars_S_const[atom2_typeid] + 
							       dockpars_qasp*fabs(KerConst->atom_charges_const[atom2_id]))*KerConst->dspars_V_const[atom1_typeid]) *
					                       dockpars_coeff_desolv*half_exp(-distance_leo*half_divide(distance_leo,25.92f));
#else	// Full precision // from 148 to 151 sec
			partial_energies[get_local_id(0)] += ((KerConst->dspars_S_const[atom1_typeid] + 
							       dockpars_qasp*fabs(KerConst->atom_charges_const[atom1_id]))*KerConst->dspars_V_const[atom2_typeid] +
					                      (KerConst->dspars_S_const[atom2_typeid] + 
							       dockpars_qasp*fabs(KerConst->atom_charges_const[atom2_id]))*KerConst->dspars_V_const[atom1_typeid]) *
					                       dockpars_coeff_desolv*exp(-distance_leo*distance_leo/25.92f);
#endif

		} // End if-statement on distance_leo 



	} // End contributor_counter for-loop

	barrier(CLK_LOCAL_MEM_FENCE);

	if (get_local_id(0) == 0)
	{
		*energy = partial_energies[0];

		for (contributor_counter=1; 
		     contributor_counter<NUM_OF_THREADS_PER_BLOCK;
		     contributor_counter++)
		{
			*energy += partial_energies[contributor_counter];
		}
	}


}
