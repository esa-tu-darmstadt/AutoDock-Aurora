#include "auxiliary.h"

void dummy_ada(
			ushort				DockConst_max_num_of_iters,
			float				DockConst_rho_lower_bound,
			uchar				DockConst_num_of_genes,
			float               DockConst_gridsize_x_minus1,
			uint				pop_size
)
{
	printf("\n");
	printf("Starting <dummy ada> ... \n");
	printf("\n");
	printf("LS: DockConst_max_num_of_iters: %u\n",		DockConst_max_num_of_iters);
	printf("LS: DockConst_rho_lower_bound: %f\n",      	DockConst_rho_lower_bound);
	printf("LS: DockConst_num_of_genes: %u\n",  	   	DockConst_num_of_genes);
	printf("LS: DockConst_gridsize_x_minus1: %f\n",  	DockConst_gridsize_x_minus1);
	printf("LS: pop_size: %u\n",           				pop_size);
}