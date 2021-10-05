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

#ifndef PROCESSRESULT_H_
#define PROCESSRESULT_H_

#include "ext_headers.h"
#include "processligand.h"
#include "getparameters.h"

#define PRINT1000(file, x) fprintf(file,  ((fabs((x)) >= 0.0) && ((fabs(x)) <= 1000.)) ? "%+7.2f" : "%+11.2e" , (x));

typedef struct
{
	Liganddata reslig_realcoord;
	float 	   interE;
	float 	   interE_elec;
	float      intraE;
	float      peratom_vdw  [MAX_NUM_OF_ATOMS];
	float      peratom_elec [MAX_NUM_OF_ATOMS];
	float      rmsd_from_ref;
	float      rmsd_from_cluscent;
	int        clus_id;
	int        clus_subrank;
	int        run_number;
} Ligandresult;


void arrange_result(    float*    final_population,
		        float*    energies,
		    const int 	  pop_size);

void write_basic_info(            FILE* fp,
		      const Liganddata* ligand_ref,
		      const Dockpars*   mypars,
		      const Gridinfo*   mygrid,
		      const int*        argc,
		           char**       argv);

void write_basic_info_dlg(	      FILE* fp,
			  const Liganddata* ligand_ref,
			  const Dockpars*   mypars,
			  const Gridinfo*   mygrid,
			  const int*        argc,
			  	char**      argv);

void make_resfiles(	      
				float* 		final_population,
				float* 		energies,
		const 	Liganddata* ligand_ref,
        const 	Liganddata* ligand_from_pdb,
		const 	Liganddata* ligand_xray,
		const 	Dockpars*   mypars,
				int  		evals_performed,
            	int 		generations_used,
        const 	Gridinfo*   mygrid,
        const 	float*      grids,
        		float* 		cpu_ref_ori_angles,
        const 	int* 	    argc,
				char** 		argv,
		   		int  		debug,
		   		int  		run_cnt,
		   		Ligandresult* best_result
					);

void cluster_analysis(     Ligandresult myresults [],
		                    int num_of_runs,
		                  char* report_file_name,
		      const Liganddata* ligand_ref,
		      const Dockpars* mypars,
		      const Gridinfo* mygrid,
		      const int*      argc,
		            char**    argv,
		      const double    docking_avg_runtime,
		      const double    program_runtime);

void clusanal_gendlg(Ligandresult myresults [],
		                  int  num_of_runs,
		     const Liganddata* ligand_ref,
		     const Dockpars*   mypars,
                     const Gridinfo*   mygrid,
		     const int*        argc,
                               char**  argv,
                     const double docking_avg_runtime,
		     const double program_runtime);

#endif /* PROCESSRESULT_H_ */
