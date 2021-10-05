/*

AutoDock-Aurora, a vectorized implementation of AutoDock 4.2 for the NEC SX-Aurora TSUBASA
Copyright (C) 2021 TU Darmstadt, Embedded Systems and Applications Group, Germany. All rights reserved.
For some of the code, Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

AutoDock is a Trade Mark of the Scripps Research Institute.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <stdio.h>
#include <stdlib.h>
#include "rand_gen.h"

#if defined (ENABLE_TRACE)
#include <ftrace.h>
#endif

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
  ftrace_region_begin("randf_vec");
#endif
#pragma omp critical
  {
    asl_random_generate_s(rng, num, val);
  }
#if defined (ENABLE_TRACE)
  ftrace_region_end("randf_vec");
#endif
}
