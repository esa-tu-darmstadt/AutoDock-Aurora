#include "auxiliary.c"

// --------------------------------------------------------------------------
// Calculation of the pose of the ligand according to the input genotype.
// --------------------------------------------------------------------------
void calc_pc (
	const	int*	restrict	PC_rotlist,
	const	int* 	restrict	PC_subrotlist_1,
	const	int* 	restrict	PC_subrotlist_2,
	const	int* 	restrict	PC_subrotlist_3,
	const	int* 	restrict	PC_subrotlist_4,
	const	int* 	restrict	PC_subrotlist_5,
	const	int* 	restrict	PC_subrotlist_6,
	const	int* 	restrict	PC_subrotlist_7,
	const	float*	restrict	PC_ref_coords_x,
	const	float*	restrict	PC_ref_coords_y,
	const	float*	restrict	PC_ref_coords_z,
	const	float*	restrict	PC_rotbonds_moving_vectors,
	const	float*	restrict	PC_rotbonds_unit_vectors,
	const	float*	restrict	PC_ref_orientation_quats,
			uint				DockConst_rotbondlist_length,
			uchar				DockConst_num_of_genes,
			uint				Host_RunId,

	const	float*	restrict	genotype,
			float*	restrict 	local_coords_x,
			float*	restrict	local_coords_y,
			float*	restrict	local_coords_z
)
{
	#if defined (PRINT_ALL_PC) 
	printf("\n");
	printf("Starting <pose calculation> ... \n");
	printf("\n");
	printf("\t%-40s %u\n", "DockConst_rotbondlist_length: ",	DockConst_rotbondlist_length);
	printf("\t%-40s %u\n", "DockConst_num_of_genes: ",        	DockConst_num_of_genes);
	printf("\t%-40s %u\n", "Host_RunId: ",        				Host_RunId);
	#endif

	#if defined (ENABLE_TRACE)
	ftrace_region_begin("PC_GENOTYPES_LOOP");
	#endif

	float local_genotype [ACTUAL_GENOTYPE_LENGTH];

	for (uint i = 3; i < DockConst_num_of_genes; i++) {
		local_genotype [i] = genotype[i] * DEG_TO_RAD;
	}

	#if defined (ENABLE_TRACE)
	ftrace_region_end("PC_GENOTYPES_LOOP");
	#endif

	float phi = local_genotype[3];
	float theta = local_genotype[4];
	float genrotangle = local_genotype[5];

	float sin_theta = sin(theta);
	float cos_theta = cos(theta);
	float genrot_unitvec[3];
	genrot_unitvec[0] = sin_theta*cos(phi);
	genrot_unitvec[1] = sin_theta*sin(phi);
	genrot_unitvec[2] = cos_theta;

	#if defined (ENABLE_TRACE)
	ftrace_region_begin("PC_MAIN_LOOP");
	#endif

	calc_pc_main_loop (
    	PC_rotlist,
    	PC_ref_coords_x,
		PC_ref_coords_y,
		PC_ref_coords_z,
		PC_rotbonds_moving_vectors,
		PC_rotbonds_unit_vectors,
		PC_ref_orientation_quats,
        DockConst_rotbondlist_length,
        Host_RunId,
    	genotype,
    	local_coords_x,
		local_coords_y,
		local_coords_z,
		genrot_unitvec,
        local_genotype
	);

	#if defined (ENABLE_TRACE)
	ftrace_region_end("PC_MAIN_LOOP");
	#endif

	#if defined (PRINT_ALL_PC) 
	printf("\n");
	printf("Finishing <pose calculation>\n");
	printf("\n");
	#endif
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
