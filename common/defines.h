#ifndef DEFINES_H_
#define DEFINES_H_

// TODO: check with current ones in AutoDock-GPU

// -----------------------------------------------------------------------------
// DOCKING CONSTANTS
// -----------------------------------------------------------------------------
#define RLIST_ATOMID_MASK    0x000000FF
#define RLIST_RBONDID_MASK   0x0000FF00
#define RLIST_RBONDID_SHIFT  8
#define RLIST_FIRSTROT_MASK  0x00010000
#define RLIST_GENROT_MASK    0x00020000
#define RLIST_DUMMY_MASK     0x00040000
//#define DEG_TO_RAD 	     0.0174533f
#define DEG_TO_RAD	     0.0174532925f

//#define COEFF_VDW 0.1662f
//#define COEFF_HB 0.1209f
//#define COEFF_ELEC 46.6881f
//#define COEFF_DESOLV 0.1322f
//#define QASP_AD 0.01097f
//#define QASP_DS 0.00679f

#define ATYPE_NUM				22
#define MAX_NUM_OF_ATOMS		256
#define MAX_NUM_OF_ATYPES		14
#define MAX_INTRAE_CONTRIBUTORS	8128
#define MAX_NUM_OF_ROTATIONS	4096
#define MAX_NUM_OF_ROTBONDS		32
#define MAX_NUM_GENES			(MAX_NUM_OF_ROTBONDS + 6)

//TODO: DECIDE WHICH IS BETTER
#define MAX_POPSIZE 			2048
//#define MAX_POPSIZE 			256
#define MAX_NUM_OF_RUNS			1000

// Must be bigger than MAX_NUM_OF_ROTBONDS+6
#define GENOTYPE_LENGTH_IN_GLOBMEM 64
#define ACTUAL_GENOTYPE_LENGTH	(MAX_NUM_OF_ROTBONDS+6)

#define LS_EXP_FACTOR 	2.0f
#define LS_CONT_FACTOR 	0.5f

//macro for a%b where b=2^N
//#define MOD2N(a, b) (a&(b-1))

// gcc linear congruential generator constants
#define RAND_A 		1103515245u
#define RAND_C 		12345u
//WARNING: it is supposed that unsigned int is 32 bit long
#define MAX_UINT 	4294967296.0f
#define MAX_ONE_FACTOR	(0.999999f / MAX_UINT)

// Macro for trilinear interpolation
#define TRILININTERPOL(cube, weights) (cube[0][0][0]*weights[0][0][0] +cube[1][0][0]*weights[1][0][0] +	\
				       cube[0][1][0]*weights[0][1][0] +cube[1][1][0]*weights[1][1][0] + \
				       cube[0][0][1]*weights[0][0][1] +cube[1][0][1]*weights[1][0][1] + \
				       cube[0][1][1]*weights[0][1][1] +cube[1][1][1]*weights[1][1][1])

// Constants for dielectric term of the
// electrostatic component of the intramolecular energy/gradient
#define DIEL_A					-8.5525f
#define DIEL_WAT				78.4f
#define DIEL_B					(DIEL_WAT - DIEL_A)	// 86.9525f
#define DIEL_LAMBDA				0.003627f
#define DIEL_H					DIEL_LAMBDA
#define DIEL_K					7.7839f
#define DIEL_B_TIMES_H			(DIEL_B * DIEL_H)	// 0.315376718f
#define DIEL_B_TIMES_H_TIMES_K	(DIEL_B_TIMES_H * DIEL_K)	// 2.454860831f

// Constant for gradients
#define PI_FLOAT						3.14159265f
#define PI_TIMES_2						(2.0f * PI_FLOAT)

#define INFINITESIMAL_RADIAN			(1e-3f)
#define HALF_INFINITESIMAL_RADIAN		(0.5f * INFINITESIMAL_RADIAN)
#define INV_INFINITESIMAL_RADIAN		(1.0f / INFINITESIMAL_RADIAN)
#define COS_HALF_INFINITESIMAL_RADIAN 	cosf(HALF_INFINITESIMAL_RADIAN)
#define SIN_HALF_INFINITESIMAL_RADIAN 	sinf(HALF_INFINITESIMAL_RADIAN)
#define inv_angle_delta					(500.0f / PI_FLOAT)

#endif /* DEFINES_H_ */
