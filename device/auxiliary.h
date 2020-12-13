#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include <stdio.h>
#include <stdint.h>

#include "defines.h"

#include "math.h"

#if defined (ENABLE_TRACE)
#include <ftrace.h>
#endif

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
static inline
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

static inline
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
 * Pose calculation 
 * -----------------------------------------------
 * */

static inline
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

static inline
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

static inline
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

static inline
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

static inline
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
static inline
float esa_sqrt(const float x){
	//uint i = as_uint(x);
	unsigned int i = *(unsigned int*) &x;    	
	i  += 127 << 23;		// adjust bias   	
	i >>= 1; 				// approximation of square root 	
	return *(float*) &i; 	// return as_float(i);
}  

#endif /* AUXILIARY_H_ */
