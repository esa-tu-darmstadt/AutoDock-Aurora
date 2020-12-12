#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include <stdio.h>
#include <stdint.h>

#include "defines.h"

#include "math.h"

#if defined (ENABLE_TRACE)
#include <ftrace.h>
#endif

// should build this as separate object...
#include "rand_gen.c"

/* 
 * -----------------------------------------------
 * Overall
 * -----------------------------------------------
 * */
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int  uint;

typedef enum {False, True} boolean;

/*
 * -----------------------------------------------
 * Map the argument into the interval
 * 0 - 180, or 0 - 360 in sexagesimal
 * -----------------------------------------------
 * */
float map_angle_180 (float angle)
{
	float x = angle;
/*
	while (x < 0.0f)   { x += 180.0f; }
	while (x > 180.0f) { x -= 180.0f; }
	return x;
*/
	x = fmod(x, 180.0f);
	if (x < 0.0f)
		x += 180.0f;
	return x;
}

float map_angle_360 (float angle) {
	float x = angle;
/*
	while (x < 0.0f)	{ x += 360.0f; }
	while (x > 360.0f)	{ x -= 360.0f; }
	return x;
*/
	x = fmod(x, 360.0f);
	if (x < 0.0f)
		x += 360.0f;
	return x;
}

/*
 * -----------------------------------------------
 * Random
 * -----------------------------------------------
 * */

/*
 * Scalar implementation rand().
 * This kernel subfunction generates a random int
 * with a linear congruential generator (LCG).
 * It uses the gcc LCG constants.
 */
unsigned int rand(unsigned int* input) {
  unsigned int rand;

  // Calculating next state
#if defined (REPRO)
  rand = 1;
#else
  rand = (RAND_A * input[0] + RAND_C);
#endif

  // Saving next state to memory
  input[0] = rand;

  return rand;
}

/*
 * Scalar implementation randf().
 * This kernel subfunction generates a random floatq
 * greater than (or equal to) 0 and less than 1.
 * It uses rand() function.
 */
float randf(unsigned int* input) {
  float randf;

  // State will be between 0 and 1
#if defined (REPRO)
  randf = 0.55f;
#else
/*
	randf = (rand(input) / MAX_UINT) *0.999999f;
*/
	randf = rand(input) * MAX_ONE_FACTOR;
#endif

  return randf;
}

/* 
 * -----------------------------------------------
 * Pose calculation 
 * -----------------------------------------------
 * */

float esa_dot3(float a[3], float b[3]) {

	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);

/*
	float tmp[3];
	for (uint i = 0; i < 3; i++) {
		tmp[i] = a[i] * b[i];
	}
	return (tmp[0] + tmp[1] + tmp[2]);
*/
}

float esa_dot3_e(float a1, float a2, float a3, float b1, float b2, float b3) {

	return (a1*b1 + a2*b2 + a3*b3);

/*
	float tmp[3];
	for (uint i = 0; i < 3; i++) {
		tmp[i] = a[i] * b[i];
	}
	return (tmp[0] + tmp[1] + tmp[2]);
*/
}

float esa_dot4(float a[4], float b[4]) {
	
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]);

/*
	float tmp[4];
	for (uint i = 0; i < 4; i++) {
		tmp[i] = a[i] * b[i];
	}
	return (tmp[0] + tmp[1] + tmp[2] + tmp[3]);
*/
}

float esa_dot4_e(float a1, float a2, float a3, float a4, float b1, float b2, float b3, float b4) {

	return (a1*b1 + a2*b2 + a3*b3 + a4*b4);

/*
	float tmp = 0.0;
	for (uint i = 0; i < 4; i++) {
		tmp += a[i] * b[i];
	}
	return tmp;
*/
}

/* 
 * -----------------------------------------------
 * Intermolecular
 * -----------------------------------------------
 * */

int esa_ceil (float fp) {
  return ((-1) * floor(-1*fp));
}

/* 
 * -----------------------------------------------
 * Intramolecular
 * -----------------------------------------------
 * */

// sqrt7 
//https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi
float esa_sqrt(const float x){
	//uint i = as_uint(x);
	unsigned int i = *(unsigned int*) &x;    	
	i  += 127 << 23;		// adjust bias   	
	i >>= 1; 				// approximation of square root 	
	return *(float*) &i; 	// return as_float(i);
}  

#endif /* AUXILIARY_H_ */
