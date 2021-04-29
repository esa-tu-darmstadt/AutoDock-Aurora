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
