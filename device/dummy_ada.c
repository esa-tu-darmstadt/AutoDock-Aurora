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
