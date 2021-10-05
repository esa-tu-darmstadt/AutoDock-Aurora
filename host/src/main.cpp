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

#include "ext_headers.h"
#include "processgrid.h"
#include "processligand.h"
#include "getparameters.h"
#include "performdocking.h"

int main(int argc, char* argv[])
{

	//=======================================================================
	// Docking Algorithm
	//=======================================================================

	Gridinfo mygrid;
	Liganddata myligand_init;
	Dockpars mypars;
	float* 	floatgrids;
	FILE*	fp;
	char 	report_file_name [256];

	clock_t clock_start_program;
	clock_start_program = clock();

	// ------------------------
	// Time measurement
	double num_sec, num_usec, elapsed_sec;
	timeval time_start,time_end;
	gettimeofday(&time_start,NULL);
	// ------------------------

	//------------------------------------------------------------
	// Capturing names of grid parameter file and ligand pdbqt file
	//------------------------------------------------------------

	// Filling the filename and coeffs fields of mypars according to command line arguments
	if (get_filenames_and_ADcoeffs(&argc, argv, &mypars) != 0)
		return 1;

	//------------------------------------------------------------
	// Processing receptor and ligand files
	//------------------------------------------------------------

	// Filling mygrid according to the gpf file
	if (get_gridinfo(mypars.fldfile, &mygrid) != 0)
		return 1;

	//allocating CPU memory for floatgrids
	size_t size_fgrid_nelems = (mygrid.num_of_atypes+2) * mygrid.size_xyz[0] * mygrid.size_xyz[1] * mygrid.size_xyz[2];

	// Filling the atom types filed of myligand according to the grid types
	if (init_liganddata(mypars.ligandfile, &myligand_init, &mygrid) != 0)
		return 1;

	// Filling myligand according to the pdbqt file
	if (get_liganddata(mypars.ligandfile, &myligand_init, mypars.coeffs.AD4_coeff_vdW, mypars.coeffs.AD4_coeff_hb) != 0)
		return 1;

	//Reading the grid files and storing values in the memory region pointed by floatgrids
	if (get_gridvalues_f(&mygrid, &floatgrids) != 0)
		return 1;

	//------------------------------------------------------------
	// Capturing algorithm parameters (command line args)
	//------------------------------------------------------------
	get_commandpars(&argc, argv, &(mygrid.spacing), &mypars);

	Liganddata myxrayligand;
	Gridinfo   mydummygrid;
	// if -lxrayfile provided, then read xray ligand data
	if (mypars.given_xrayligandfile == true) {
			if (init_liganddata(mypars.xrayligandfile, &myxrayligand, &mydummygrid) != 0)
				return 1;

			if (get_liganddata(mypars.xrayligandfile, &myxrayligand, mypars.coeffs.AD4_coeff_vdW, mypars.coeffs.AD4_coeff_hb) != 0)
				return 1;
	}

	//------------------------------------------------------------
	// Calculating energies of reference ligand if required
	//------------------------------------------------------------
	if (mypars.reflig_en_reqired == 1) {
		print_ref_lig_energies_f(
					myligand_init,
					mypars.smooth,
					mygrid,
					floatgrids,
					mypars.coeffs.scaled_AD4_coeff_elec,
					mypars.coeffs.AD4_coeff_desolv,
					mypars.qasp);
	}

	//------------------------------------------------------------
	// Starting Docking
	//------------------------------------------------------------
	if (docking_with_aurora(&mygrid, floatgrids, &mypars, &myligand_init, &myxrayligand, &argc, argv, clock_start_program) != 0)
		return 1;

	free(floatgrids);

	// ------------------------
	// Time measurement
	gettimeofday(&time_end,NULL);
	num_sec     = time_end.tv_sec  - time_start.tv_sec;
	num_usec    = time_end.tv_usec - time_start.tv_usec;
	elapsed_sec = num_sec + (num_usec/1000000);
	printf("Program run time %.3f sec\n\n",elapsed_sec);

	// Append time information to .dlg file
	strcpy(report_file_name, mypars.resname);
	strcat(report_file_name, ".dlg");
	fp = fopen(report_file_name, "a");
	fprintf(fp, "\n\n\nProgram run time %.3f sec\n", elapsed_sec);
	fclose(fp);
	// ------------------------

	return 0;
}
