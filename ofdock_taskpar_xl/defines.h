#ifndef DEFINES_H_
#define DEFINES_H_

// -----------------------------------------------------------------------------
// DIFFERENT WORK GROUP SIZES
// -----------------------------------------------------------------------------
/*
//NUM_WORKITEMS=-DN32WI 
#if defined (N32WI) 	
	#define NUM_OF_THREADS_PER_BLOCK 32  // ofdock_amd_32wi 
#elif defined (N16WI) 	
	#define NUM_OF_THREADS_PER_BLOCK 16  // ofdock_amd_16wi 
#elif defined (N64WI) 	
	#define NUM_OF_THREADS_PER_BLOCK 64  // ofdock_amd_64wi 
#elif defined (N128WI) 	
	#define NUM_OF_THREADS_PER_BLOCK 128 // ofdock_amd_128wi 
#else 	
	#define NUM_OF_THREADS_PER_BLOCK 1  // ofdock_amd_32wi 
#endif
*/

#define NUM_OF_THREADS_PER_BLOCK 1
// -----------------------------------------------------------------------------
// DOCKING CONSTANTS
// -----------------------------------------------------------------------------
#define RLIST_ATOMID_MASK    0x000000FF
#define RLIST_RBONDID_MASK   0x0000FF00
#define RLIST_RBONDID_SHIFT  8
#define RLIST_FIRSTROT_MASK  0x00010000
#define RLIST_GENROT_MASK    0x00020000
#define RLIST_DUMMY_MASK     0x00040000
#define DEG_TO_RAD 	     0.0174533f

//#define COEFF_VDW 0.1662f
//#define COEFF_HB 0.1209f
//#define COEFF_ELEC 46.6881f
//#define COEFF_DESOLV 0.1322f
//#define QASP_AD 0.01097f
//#define QASP_DS 0.00679f

#define MAX_NUM_OF_ATOMS 	90
#define MAX_NUM_OF_ATYPES 	14
#define MAX_INTRAE_CONTRIBUTORS 8128
#define MAX_NUM_OF_ROTATIONS 	4096
#define MAX_NUM_OF_ROTBONDS 	32
//#define MAX_POPSIZE 		2048
#define MAX_POPSIZE 		256
#define MAX_NUM_OF_RUNS 	100

// Must be bigger than MAX_NUM_OF_ROTBONDS+6
#define GENOTYPE_LENGTH_IN_GLOBMEM 64
#define ACTUAL_GENOTYPE_LENGTH	(MAX_NUM_OF_ROTBONDS+6)

#define LS_EXP_FACTOR 	2.0f
#define LS_CONT_FACTOR 	0.5f

// Mask for genotypes
#define MASK_GENOTYPE 0x3F

//macro for a%b where b=2^N
//#define MOD2N(a, b) (a&(b-1))

// gcc linear congruential generator constants
#define RAND_A 		1103515245u
#define RAND_C 		12345u
//WARNING: it is supposed that unsigned int is 32 bit long
#define MAX_UINT 	4294967296.0f

// Macro for capturing grid values
// Original
#define GETGRIDVALUE(mempoi,gridsize_x,gridsize_y,gridsize_z,t,z,y,x)   *(mempoi + gridsize_x*(y + gridsize_y*(z + gridsize_z*t)) + x)

// Optimization 1
//#define GETGRIDVALUE_OPT(mempoi,gridsize_x,gridsize_y,mul_tmp,z,y,x)   *(mempoi + gridsize_x*(y + gridsize_y*(z + mul_tmp)) + x)

// Optimization 2
// It is done in the kernel code

// Macro for trilinear interpolation
#define TRILININTERPOL(cube, weights) (cube[0][0][0]*weights[0][0][0] +cube[1][0][0]*weights[1][0][0] +	\
				       cube[0][1][0]*weights[0][1][0] +cube[1][1][0]*weights[1][1][0] + \
				       cube[0][0][1]*weights[0][0][1] +cube[1][0][1]*weights[1][0][1] + \
				       cube[0][1][1]*weights[0][1][1] +cube[1][1][1]*weights[1][1][1])

// -----------------------------------------------------------------------------
// IMPROVEMENTS OVER PECHAN'S IMPLEMENTATION
// -----------------------------------------------------------------------------

//#define NATIVE_PRECISION	// if disabled, then full precision, 
				// improvement: 152 (disabled) to 75 (enabled) sec
//#define HALF_PRECISION		// if disabled, then full precision,
				//no improvement: 152 (disabled) to 153 (enabled) sec

#define ASYNC_COPY		// if disabled, then original approach of copy, 
				// no improvement: 76 sec both
#define IMPROVE_GRID		// if disabled, then original approach to calculate grid contribution
				// improvement from 78 (disabled) to 76 (enabled) sec

#define MAPPED_COPY	// if disabled, then mem copy inside main computation loop
			// doesn't use DMA engine
			// improvement: 77 (disabled) to 76 (enabled) sec 
//#define REPRO
// -----------------------------------------------------------------------------

// Enable to debug kernels using printf
// Disable them to save RAM resources
//#define DEBUG_ACTIVE_KERNEL
//#define DEBUG_KRNL_GA
//#define DEBUG_KRNL_IC
//
//
//#define DEBUG_KRNL_CONFORM
//#define DEBUG_KRNL_INTERE
//#define DEBUG_KRNL_INTRAE
//#define DEBUG_KRNL_STORE


#endif /* DEFINES_H_ */
