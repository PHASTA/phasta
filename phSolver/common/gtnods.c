#include <stdio.h>
#include <mpi.h>

#include "common_c.h"

#include <FCMangle.h>  
#define gtnods FortranCInterface_GLOBAL_(gtnods, GTNODS) 

void gtnods()
{

	if(workfc.numpe > 1) {

		int irecvcount, ierr;
		long long numvec, loc_nshgt;

                irecvcount = 1;
		numvec = (long long) newdim.nshg0;

//		printf("Local number of modes = %ld %d %d\n",numvec,newdim.nshg0, sizeof(newdim.nshg0));

		ierr = MPI_Allreduce(&numvec, &loc_nshgt, irecvcount,
                                MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		newdim.nshgt = loc_nshgt;	
	}
	else {
		newdim.nshgt = (long long) conpar.nshg;
	}

	if (workfc.myrank == workfc.master) {
		printf("Total number of modes = %ld\n",newdim.nshgt);
 	}

}
