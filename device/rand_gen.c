#include <stdio.h>
#include <asl.h>

static asl_random_t rng;

void randf_vec_init()
{
  int rc;

  asl_library_initialize();

  rc = asl_random_create(&rng, ASL_RANDOMMETHOD_MT19937_64);
  if (rc != 0)
    printf("ERROR: asl_random_create rc=%d\n", rc);

  rc = asl_random_distribute_uniform(rng);
  if (rc != 0)
    printf("ERROR: asl_random_distribute_uniform returned rc=%d\n", rc);
}

void randf_vec_fini()
{
    asl_random_destroy(rng);
    asl_library_finalize();
}

void rand_vec(float *val, int num)
{
#ifdef _OPENMP
#pragma omp critical
#endif
  asl_random_generate_s(rng, num, val);
}


#if 1
int main(int argc, char *argv[])
{
  int n = 100;
  float val[n];
  
  randf_vec_init();
  rand_vec(val, n);
  for (int i = 0; i < n; i++)
    printf("%3d  %f\n", i, val[i]);
  randf_ve_fini();
}
#endif
