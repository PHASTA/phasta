#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>
#ifdef __bgq__
#include <hwi/include/bqc/A2_inlines.h>
#endif
#include "FCMangle.h"

static uint64_t start;
static uint64_t total;


#define cycle_count_stop FortranCInterface_GLOBAL_(cycle_count_stop, CYCLE_COUNT_STOP)
#define cycle_count_start FortranCInterface_GLOBAL_(cycle_count_start, CYCLE_COUNT_START)
#define cycle_count_print FortranCInterface_GLOBAL_(cycle_count_print, CYCLE_COUNT_PRINT)

void cycle_count_stop()
{
#ifdef __bgq__
	uint64_t tb = GetTimeBase();
	if(__builtin_expect(tb < start, 0))
	{
		total += tb + (UINT64_MAX-start);
	}
	else
	{
		total += tb-start;
	}
#else
	total = 0; //FIEX ME for non BGQ env
#endif
}
void cycle_count_start()
{
#ifdef __bgq__
	start = GetTimeBase();
#else
	start = 0; //FIEX ME for non BGQ env
#endif
}

void cycle_count_print()
{
        uint64_t iresult; 
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Allreduce(&total, &iresult, 1, MPI_LONG_LONG_INT, 
			MPI_MAX, MPI_COMM_WORLD);
        if(rank == 0)
	printf("fillsparse : %llu cycles\n", iresult);
	total = 0;
}
