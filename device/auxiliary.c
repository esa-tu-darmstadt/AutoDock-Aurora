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


#endif /* AUXILIARY_H_ */