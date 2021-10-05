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

#ifndef PERFORMDOCKING_H_
#define PERFORMDOCKING_H_

#include "ext_headers.h"
#include "processgrid.h"
#include "miscellaneous.h"
#include "processligand.h"
#include "getparameters.h"
#include "calcenergy.h"	
#include "processresult.h"

#include "veo_api_wrappers.h"

#define sizeKB(x) x/1000.0f

#define ELAPSEDSECS(stop,start) ((float) stop-start)/((float) CLOCKS_PER_SEC)

int docking_with_aurora(
	const	Gridinfo*	mygrid,
	/*const*/ float*	cpu_floatgrids,
			Dockpars*	mypars,
	const	Liganddata*	myligand_init,
	const 	Liganddata* myxrayligand,
	const	int*		argc,
			char**		argv,
			clock_t		clock_start_program
);

#endif /* PERFORMDOCKING_H_ */
