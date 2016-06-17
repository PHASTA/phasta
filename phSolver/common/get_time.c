#include <stdio.h>
#include <sys/types.h>
#include <mpi.h>
#include <stdint.h>
#include <assert.h>
#ifdef __bgq__
#include <hwi/include/bqc/A2_inlines.h>
#endif
#include "FCMangle.h"

#define get_time FortranCInterface_GLOBAL_(get_time, GET_TIME)
#define get_max_time_diff FortranCInterface_GLOBAL_(get_max_time_diff, GET_MAX_TIME_DIFF)
static double multiplier;

void get_max_time_diff(uint64_t* first, uint64_t* last, uint64_t* c_first, uint64_t* c_last, char* lbl)
{
	uint64_t tmp = ((*last)-(*first));
	uint64_t iresult;
	double result;
	uint64_t c_result;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        assert(tmp > 0);
	MPI_Allreduce(&tmp, &iresult, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
	result = iresult*multiplier;
	if(__builtin_expect(*c_last < *c_first, 0))
	{
		tmp = *c_last + (UINT64_MAX - *c_first);
	}
	else
	{
		tmp = ((*c_last)-(*c_first));
	}
	MPI_Allreduce(&tmp, &iresult, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
	c_result = iresult;
	if(rank == 0)
		printf("%s: %.12g seconds, %llu cycles\n", lbl, result, c_result);
}

/*
#ifdef __APPLE__
#include <mach/time.h>
#endif
*/

/*#if defined(__bgq__) || defined(__APPLE__)*/
#if (1)
#include <sys/time.h>
void get_time(uint64_t* rv, uint64_t* cycle)
{
	struct timeval time;
	int ret;
	ret = gettimeofday(&time, NULL);
	if(ret != 0) perror("gettimeofday failed: ");
	*rv = ((time.tv_sec*1000000)+(time.tv_usec));
	multiplier=0.000001;/*10.0e-6;*/
#ifdef __bgq__
	*cycle = GetTimeBase();
#endif
}
#else
#include <time.h>
void get_time(uint64_t* rv)
{
	struct timespec time;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time);
        *rv = ((time.tv_sec*1000000000)+(time.tv_nsec));
	multiplier = 0.000000001;
}
#endif
