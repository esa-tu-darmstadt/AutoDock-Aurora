#include "auxiliary.c"
#include "calc_pc_main_loop.c"

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
	const	int* 	restrict	PC_subrotlist_8,
	const	int* 	restrict	PC_subrotlist_9,
	const	int* 	restrict	PC_subrotlist_10,
	const	int* 	restrict	PC_subrotlist_11,
	const	float*	restrict	PC_ref_coords_x,
	const	float*	restrict	PC_ref_coords_y,
	const	float*	restrict	PC_ref_coords_z,
	const	float*	restrict	PC_rotbonds_moving_vectors,
	const	float*	restrict	PC_rotbonds_unit_vectors,
	const	float*	restrict	PC_ref_orientation_quats,
			uint				DockConst_rotbondlist_length,
			uint				subrotlist_1_length,
			uint				subrotlist_2_length,
			uint				subrotlist_3_length,
			uint				subrotlist_4_length,
			uint				subrotlist_5_length,
			uint				subrotlist_6_length,
			uint				subrotlist_7_length,
			uint				subrotlist_8_length,
			uint				subrotlist_9_length,
			uint				subrotlist_10_length,
			uint				subrotlist_11_length,
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

	float sin_theta = sin(theta);
	float cos_theta = cos(theta);
	float genrot_unitvec[3];
	genrot_unitvec[0] = sin_theta*cos(phi);
	genrot_unitvec[1] = sin_theta*sin(phi);
	genrot_unitvec[2] = cos_theta;

	#if defined (ENABLE_TRACE)
	ftrace_region_begin("PC_MAIN_LOOP");
	#endif

	// This was the original call processing a single large rotlist
/*
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
*/
	// The alternative version consists of splitting
	// the single rotlist into many sub-rotlists, where
	// applying ivdep is correct.
	// Note loops within sub-rotlists can be run in parallel,
	// while sub-rotlists must be processed one after another.

	if (subrotlist_1_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_1,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_1_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_2_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_2,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_2_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_3_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_3,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_3_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_4_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_4,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_4_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_5_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_5,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_5_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_6_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_6,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_6_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_7_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_7,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_7_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_8_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_8,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_8_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_9_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_9,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_9_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_10_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_10,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_10_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

	if (subrotlist_11_length > 0) {
		calc_pc_main_loop (
			PC_subrotlist_11,
			PC_ref_coords_x,
			PC_ref_coords_y,
			PC_ref_coords_z,
			PC_rotbonds_moving_vectors,
			PC_rotbonds_unit_vectors,
			PC_ref_orientation_quats,
			subrotlist_11_length,
			Host_RunId,
			genotype,
			local_coords_x,
			local_coords_y,
			local_coords_z,
			genrot_unitvec,
			local_genotype
		);
	}

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
