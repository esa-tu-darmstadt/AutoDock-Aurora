#ifndef _RAND_GEN_H_
#define _RAND_GEN_H_

#include <asl.h>

void randf_vec_init(unsigned int *seed, int seed_length);
void randf_vec_fini();
void randf_vec(float *val, int num);

#endif
