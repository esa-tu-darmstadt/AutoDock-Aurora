#ifndef AUXILIARY_H_
#define AUXILIARY_H_

/* 
 * -----------------------------------------------
 * Pose calculation 
 * -----------------------------------------------
 * */

// TODO: can be vectorized?
float esa_dot3(float a[3], float b[3]) {
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

// TODO: can be vectorized?
float esa_dot4(float a[4], float b[4]) {
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]);
}

/* 
 * -----------------------------------------------
 * Intermolecular
 * -----------------------------------------------
 * */

/* 
 * -----------------------------------------------
 * Intramolecular
 * -----------------------------------------------
 * */

// sqrt7 
//https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi
float sqrt_custom(const float x) 
{ 	//uint i = as_uint(x);	
	unsigned int i = *(unsigned int*) &x;    	
	i  += 127 << 23;		// adjust bias   	
	i >>= 1; 				// approximation of square root 	
	return *(float*) &i; 	// return as_float(i);
}  



#endif /* AUXILIARY_H_ */
