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

#include "processgrid.h"

int get_gridinfo(const char* fldfilename, Gridinfo* mygrid)
{
	FILE*  fp;
	char   tempstr [256];
	int    gpoints_even[3];
	int    recnamelen;
	double center[3];

	// ----------------------------------------------------
	// Getting full path fo the grid file
	// Getting father directory name
	//char* dir = dirname(ts1);
	//char* filename = basename(ts1);

	char* ts1 = strdup(fldfilename);
	mygrid->grid_file_path = dirname(ts1);
	// ----------------------------------------------------

	//Processing fld file
	fp = fopen(fldfilename, "r");
	if (fp == NULL)
	{
		printf("Error: can't open fld file %s!\n", fldfilename);
		return 1;
	}

	while (fscanf(fp, "%s", tempstr) != EOF)
	{
		// -----------------------------------
		// Reorder according to file *.maps.fld
		// -----------------------------------
		//Grid spacing
		if (strcmp(tempstr, "#SPACING") == 0)
		{
			fscanf(fp, "%lf", &(mygrid->spacing));
			if (mygrid->spacing > 1)
			{
				printf("Error: grid spacing is too big!\n");
				return 1;
			}
		}

		//capturing number of grid points
		if (strcmp(tempstr, "#NELEMENTS") == 0)
		{
			fscanf(fp, "%d%d%d", &(gpoints_even[0]), &(gpoints_even[1]), &(gpoints_even[2]));
			//plus one gridpoint in each dimension
			mygrid->size_xyz[0] = gpoints_even[0] + 1;
			mygrid->size_xyz[1] = gpoints_even[1] + 1;
			mygrid->size_xyz[2] = gpoints_even[2] + 1;

			//If the grid is too big, send message and change the value of truncated_size_xyz
			if ((mygrid->size_xyz [0] > 128) || (mygrid->size_xyz [1] > 128) || (mygrid->size_xyz [2] > 128))
			{
				printf("Error: each dimension of the grid must be below 128.\n");
				return 1;
			}
		}

		//Capturing center
		if (strcmp(tempstr, "#CENTER") == 0)
		{
			fscanf(fp, "%lf%lf%lf", &(center[0]), &(center[1]), &(center[2]));
		}

		//Name of the receptor and corresponding files
		if (strcmp(tempstr, "#MACROMOLECULE") == 0)
		{
			fscanf(fp, "%s", tempstr);
			recnamelen = strcspn(tempstr,".");
			tempstr[recnamelen] = '\0';
			strcpy(mygrid->receptor_name, tempstr);
		}

		// -----------------------------------
		// MISSING: similar section corresponding to
		// #GRID_PARAMETER_FILE
		// -----------------------------------
	}

	//calculating grid size
	mygrid->size_xyz_angstr[0] = (mygrid->size_xyz[0]-1)*(mygrid->spacing);
	mygrid->size_xyz_angstr[1] = (mygrid->size_xyz[1]-1)*(mygrid->spacing);
	mygrid->size_xyz_angstr[2] = (mygrid->size_xyz[2]-1)*(mygrid->spacing);

	//calculating coordinates of origo
	mygrid->origo_real_xyz[0] = center[0] - (((double) gpoints_even[0])*0.5*(mygrid->spacing));
	mygrid->origo_real_xyz[1] = center[1] - (((double) gpoints_even[1])*0.5*(mygrid->spacing));
	mygrid->origo_real_xyz[2] = center[2] - (((double) gpoints_even[2])*0.5*(mygrid->spacing));

	fclose(fp);

	return 0;
}

int get_gridvalues_f(const Gridinfo* mygrid, float** fgrids)
//The function reads the grid point values from the .map files
//that correspond to the receptor given by the first parameter.
//It allocates the proper amount of memory and stores the data there,
//which can be accessed with the fgrids pointer.
//If there are any errors, it returns 1, otherwise
//the return value is 0.
{
	int t, x, y, z;
	FILE* fp;
	char tempstr [128];
	float* mypoi;

	*fgrids = (float*) malloc((sizeof(float))*(mygrid->num_of_atypes+2)*
						  (mygrid->size_xyz[0])*
					          (mygrid->size_xyz[1])*
						  (mygrid->size_xyz[2]));
	if (*fgrids == NULL)
	{
		printf("Error: not enough memory!\n");
		return 1;
	}

	mypoi = *fgrids;

	for (t=0; t < mygrid->num_of_atypes+2; t++)
	{
		//opening corresponding .map file
		//-------------------------------------
		// Added the complete path of associated grid files.
		strcpy(tempstr,mygrid->grid_file_path);
		strcat(tempstr, "/");
		strcat(tempstr, mygrid->receptor_name);

		//strcpy(tempstr, mygrid->receptor_name);
		//-------------------------------------
		strcat(tempstr, ".");
		strcat(tempstr, mygrid->grid_types[t]);
		strcat(tempstr, ".map");
		fp = fopen(tempstr, "rb"); // fp = fopen(tempstr, "r");
		if (fp == NULL)
		{
			printf("Error: can't open %s!\n", tempstr);
			return 1;
		}

		//seeking to first data
		do    fscanf(fp, "%s", tempstr);
		while (strcmp(tempstr, "CENTER") != 0);
		fscanf(fp, "%s", tempstr);
		fscanf(fp, "%s", tempstr);
		fscanf(fp, "%s", tempstr);

		//reading values
		for (z=0; z < mygrid->size_xyz[2]; z++)
			for (y=0; y < mygrid->size_xyz[1]; y++)
				for (x=0; x < mygrid->size_xyz[0]; x++)
				{
					fscanf(fp, "%f", mypoi);
					mypoi++;
				}
	}

	return 0;
}
