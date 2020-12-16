#include <stdio.h>
#include <stdlib.h>
#include "rand_gen.h"

static asl_random_t rng;

void randf_vec_init(unsigned int *seed, int seed_length)
{
  int rc;

  asl_library_initialize();

#pragma omp critical
  {
    rc = asl_random_create(&rng, ASL_RANDOMMETHOD_MT19937_64);
    if (rc != 0)
      printf("ERROR: asl_random_create rc=%d\n", rc);

    rc = asl_random_distribute_uniform(rng);
    if (rc != 0)
      printf("ERROR: asl_random_distribute_uniform returned rc=%d\n", rc);

    rc = asl_random_initialize(rng, seed_length, seed);
    if (rc != 0)
      printf("ERROR: asl_random_initialize returned rc=%d\n", rc);
  }
}

void randf_vec_fini()
{
#pragma omp critical
  {
    asl_random_destroy(rng);
    asl_library_finalize();
  }
}

void randf_vec(float *val, int num)
{
#if defined (ENABLE_TRACE)
  ftrace_region_begin("randf_vec")
#endif
#pragma omp critical
  {
    asl_random_generate_s(rng, num, val);
  }
#if defined (ENABLE_TRACE)
  ftrace_region_end("randf_vec")
#endif
}
