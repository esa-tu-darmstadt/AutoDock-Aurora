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

#ifndef PROCESSGRID_H_
#define PROCESSGRID_H_

#include "ext_headers.h"

// libgen.h contains basename() and dir() from a fullpath name
// Specific: to open correctly grid map field fiels and associated files
// Not only during real execution but also SDAccel cpu-, hw-emulation
// http://ask.systutorials.com/681/get-the-directory-path-and-file-name-from-absolute-path-linux

#include "miscellaneous.h"

// Added for allocating fgrids
using std::vector;

#define getvalue_4Darr(mempoi, grinf, t, z, y, x)                  *(mempoi + (grinf).size_xyz[0] * (y + (grinf).size_xyz[1] * (z + (grinf).size_xyz[2]*t)) + x)
#define getvalue_4Darr_withsize(mempoi, gridsize_xyz, t, z, y, x)  *(mempoi + gridsize_xyz[0]*(y + gridsize_xyz[1] * (z + gridsize_xyz[2]*t)) + x)
//The macro helps to access the grid point values
//which were read from external files with get_gridvalues function.
//The first parameter is a pointer which points to the memory area storing the data.
//The second one is the corresponding grid info (parameter of get_gridinfo function).
//The other parameters are the type index, z, y and x coordinates of the grid point.

typedef struct
//Struct which can contain all the important informations which derives from .gpf and .xyz files.
{
	char*  grid_file_path;	  // Added to store the full path of the grid file
	char   receptor_name [64];
	int    size_xyz [3];
	double spacing;
	double size_xyz_angstr [3];
	char   grid_types [16][3];
	int    num_of_atypes;
	double origo_real_xyz [3];
} Gridinfo;

int get_gridinfo(const char*, Gridinfo*);

int get_gridvalues_f(const Gridinfo* mygrid,
						   float** fgrids);


#endif /* PROCESSGRID_H_ */
