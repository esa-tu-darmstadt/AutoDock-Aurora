#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include <stdio.h>
#include <stdint.h>

#include "defines.h"

#include "math.h"

#if defined (ENABLE_TRACE)
#include <ftrace.h>
#endif

// -----------------------------------------------
// Overall
// -----------------------------------------------
#ifndef SHORT_TYPE_NAMES_
#define SHORT_TYPE_NAMES_
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int  uint;
#endif

typedef enum {False, True} boolean;

// -----------------------------------------------
// Map the argument into the interval
// 0 - 180, or 0 - 360 in sexagesimal
// -----------------------------------------------
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

// -----------------------------------------------
// Pose calculation
// -----------------------------------------------
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

#define esa_dot3_e(a1,a2,a3,b1,b2,b3) (a1*b1 + a2*b2 + a3*b3)

static inline
float esa_dot3_e_(float a1, float a2, float a3, float b1, float b2, float b3) {

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

#define esa_dot4_e(a1,a2,a3,a4,b1,b2,b3,b4) (a1*b1 + a2*b2 + a3*b3 + a4*b4)

static inline
float esa_dot4_e_(float a1, float a2, float a3, float a4, float b1, float b2, float b3, float b4) {

	return (a1*b1 + a2*b2 + a3*b3 + a4*b4);

/*
	float tmp = 0.0;
	for (uint i = 0; i < 4; i++) {
		tmp += a[i] * b[i];
	}
	return tmp;
*/
}

// -----------------------------------------------
// Intermolecular
// -----------------------------------------------

/* https://www.codeproject.com/Tips/700780/Fast-floor-ceiling-functions */
static inline
int esa_ceil (float fp) {
#ifdef __clang__
  return ceilf(fp);
#else
  return (-floorf(-fp));
#endif
}

// -----------------------------------------------
// Intramolecular
// -----------------------------------------------

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

static inline
float esa_fabs(const float x){
	//return x >= 0.0f ? x : -x;
#ifdef __clang__
	return fabsf(x);
#else
	return fabs(x);  // is fastest!
#endif
	// https://stackoverflow.com/questions/23474796/is-there-a-fast-fabsf-replacement-for-float-in-c
	//*(unsigned int *)(&x) &= 0x7fffffff; return x;
}


static inline int __float_as_int(const float x)
{
	return *(int *)&x;
}

static inline float __int_as_float(const int x)
{
	return *(float *)&x;
}

// https://stackoverflow.com/questions/47025373/fastest-implementation-of-the-natural-exponential-function-using-sse
// max. rel. error = 3.55959567e-2 on [-87.33654, 88.72283]
static inline float esa_expf0(const float x)
{
	float a = 12102203.0f;
	int32_t b = 127 * (1 << 23) - 298765;
	return __int_as_float(__float_as_int(a * x) + b);
}

/*
static float fastExp3(const float x)  // cubic spline approximation
{
	int32_t i = (int32_t)(12102203.0f*x) + 127*(1 << 23);
	int32_t m = (i >> 7) & 0xFFFF;  // copy mantissa
	// empirical values for small maximum relative error (8.34e-5):
	i += ((((((((1277*m) >> 14) + 14825)*m) >> 14) - 79749)*m) >> 11) - 626;
	return __int_as_float(i);
}

static float fastExp4(const float x)  // quartic spline approximation
{
	int32_t i = (int32_t)(12102203.0f*x) + 127*(1 << 23);
	int32_t m = (i >> 7) & 0xFFFF;  // copy mantissa
	// empirical values for small maximum relative error (1.21e-5):
	i += (((((((((((3537*m) >> 16) + 13668)*m) >> 18) + 15817)*m) >> 14) - 80470)*m) >> 11);
	return __int_as_float(i);
}
*/

// https://forums.developer.nvidia.com/t/a-more-accurate-performance-competitive-implementation-of-expf/47528
static inline
float esa_expf(const float a){
  float f, r, j, s, t;
  int i, ia;

    // exp(a) = 2**i * exp(f); i = rintf (a / log(2))
  j = (1.442695f * a + 12582912.f) - 12582912.f; // 0x1.715476p0, 0x1.8p23
  f = -6.93145752e-1f * j + a;                   // -0x1.62e400p-1  // log_2_hi 
  f = -1.42860677e-6f * j + f;                   // -0x1.7f7d1cp-20 // log_2_lo 
  i = (int)j;
  // approximate r = exp(f) on interval [-log(2)/2, +log(2)/2]
  r =            1.37805939e-3f;  // 0x1.694000p-10
  r = r * f + 8.37312452e-3f; // 0x1.125edcp-7
  r = r * f + 4.16695364e-2f; // 0x1.555b5ap-5
  r = r * f + 1.66664720e-1f; // 0x1.555450p-3
  r = r * f + 4.99999851e-1f; // 0x1.fffff6p-2
  r = r * f + 1.00000000e+0f; // 0x1.000000p+0
  r = r * f + 1.00000000e+0f; // 0x1.000000p+0
  // exp(a) = 2**i * r;
  ia = (i > 0) ?  0 : 0x83000000;
  s = __int_as_float (0x7f000000 + ia);
  t = __int_as_float ((i << 23) - ia);
  r = r * s;
  r = r * t;
  // handle special cases: severe overflow / underflow
  if (fabsf (a) >= 104.0f) r = __int_as_float ((__float_as_int (a) > 0) ? 0x7f800000 : 0);
  return r;
}

// -----------------------------------------------
// Gradients
// -----------------------------------------------

static inline
void esa_cross3_e_(float a1, float a2, float a3, float b1, float b2, float b3, float* c1, float* c2, float* c3) {
	*c1 = (a2 * b3) - (a3 * b2);
	*c2 = (a3 * b1) - (a1 * b3);
	*c3 = (a1 * b2) - (a2 * b1);
}

#define esa_length3_e(a1,a2,a3) esa_sqrt(a1*a1 + a2*a2 + a3*a3)

static inline
void esa_normalize3_e_(float a1, float a2, float a3, float* b1, float* b2, float* b3) {
	float inv_len = 1.0f / esa_length3_e(a1, a2, a3);
	*b1 = a1 * inv_len;
	*b2 = a2 * inv_len;
	*b3 = a3 * inv_len;
}

#endif /* AUXILIARY_H_ */
