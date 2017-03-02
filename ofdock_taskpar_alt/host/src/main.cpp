/*
 * (C) 2013. Evopro Innovation Kft.
 *
 * main.c
 *
 *  Created on: 2008.09.04.
 *      Author: pechan.imre
 */

#include <stdio.h>
#include <stdlib.h>


#include "processgrid.h"
#include "miscallenous.h"
#include "processresult.h"
#include "processligand.h"
#include "getparameters.h"
#include "performdocking.h"

int main(int argc, char* argv [])
{

	Gridinfo mygrid;
	Liganddata myligand_init;
	Dockpars mypars;
	
/*
	double* floatgrids;
*/
	float* floatgrids;

	double start_program, stop_program;
	double dock_exectime_sum;
	Ligandresult* result_ligands;
	//char report_file_name [30];


	start_program = timer_gets();

	//Capturing names of grid parameter file and ligand pdbqt file

	if (get_filenames_and_ADcoeffs(&argc, argv, &mypars) != 0)		//Filling the filename and coeffs fields of mypars according to command line arguments
		return 1;

	//Processing receptor and ligand files

	if (get_gridinfo(mypars.fldfile, &mygrid) != 0)		//Filling mygrid according to the gpf file
		return 1;

	if (init_liganddata(mypars.ligandfile, &myligand_init, &mygrid) != 0)		//Filling the atom types filed of myligand according to the grid types
		return 1;

	if (get_liganddata(mypars.ligandfile, &myligand_init, mypars.coeffs.AD4_coeff_vdW, mypars.coeffs.AD4_coeff_hb) != 0)	//Filling myligand according to the pdbqt file
		return 1;

	if (get_gridvalues(&mygrid, &floatgrids) != 0)		//Reading the grid files and storing values in the memory region pointed by floatgrids
		return 1;


	//Capturing algorithm parameters

	get_recligpars(&myligand_init, &mygrid, &mypars);		//Capturing required parameters of the ligand and the grid

	get_commandpars(&argc, argv, &(mygrid.spacing), &mypars);		//capturing command line arguments

	//Allocating memory for results

	result_ligands = (Ligandresult*) malloc(sizeof(Ligandresult)*(mypars.num_of_runs));
	if (result_ligands == NULL)
	{
		printf("Error: not enough memory!\n");
		return 1;
	}

	//Calculating the energies of reference ligand if required
	if (mypars.reflig_en_reqired == 1)
		print_ref_lig_energies(myligand_init, mygrid, floatgrids, mypars.coeffs.scaled_AD4_coeff_elec, mypars.coeffs.AD4_coeff_desolv, mypars.qasp);

	//Perform docking

	docking_with_cpu(&mygrid, floatgrids, &mypars, &myligand_init, &argc, argv, &dock_exectime_sum, result_ligands);

	//Perform ranked cluster analysis

	/*if (mypars.use_cpu == 1)
		strcpy(report_file_name, "clusters_cpu.txt");
	else
		strcpy(report_file_name, "clusters_fpga.txt");*/

	//cluster_analysis(result_ligands, mypars.num_of_runs, report_file_name, &myligand_init, &mypars,
	//				 &mygrid, &argc, argv, dock_exectime_sum/mypars.num_of_runs, stop_program_before_clustering-start_program);

	clusanal_gendlg(result_ligands, mypars.num_of_runs, &myligand_init, &mypars,
					&mygrid, &argc, argv);

/*
	free(floatgrids);
	free(result_ligands);
*/
	//Stopping program

	stop_program = timer_gets();

	printf("\n\nAverage run time of one run: %.3f sec\n", dock_exectime_sum/mypars.num_of_runs);
	printf("Program run time: %.3f sec\n", (stop_program-start_program));

	return 0;
}

