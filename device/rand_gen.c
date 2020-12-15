#include <stdio.h>
#include <stdlib.h>
#include "rand_gen.h"

static int nprints = 0;
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

#if 0
    printf("randf_initialize seed: ");
    for (int i=0; i<seed_length; i++)
      printf("%10u ", seed[i]);
    printf("\n\n");
#endif

    rc = asl_random_initialize(rng, seed_length, seed);
    if (rc != 0)
      printf("ERROR: asl_random_initialize returned rc=%d\n", rc);

    nprints = 0;
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
#pragma omp critical
  {
    asl_random_generate_s(rng, num, val);
  }
#if 0
  if (nprints++ < 10) {
    printf("randf_vec: ");
    for (int i=0; i<10; i++)
      printf("%12.5e ", val[i]);
    printf("\n\n");
  }
#endif
}

