#ifdef __bgq__
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpix.h>

void print_bgq_net_info()
{
	int rank,i,iscomplete;
	MPIX_Hardware_t hw;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);
	MPIX_Hardware(&hw);
	unsigned int cores = 1;
	if(rank == 0)
	{
		//printf("Detected %u physical cores\n", hw.psize);
		//FIXME: above broken when oversubscribed
		for(i=0;i<5;i++)
		{
		       if(hw.Size[i] > 0)	
				cores = cores*hw.Size[i];
		}
		cores = cores*16;
		printf("Detected %u physical cores\n", cores);
		printf("%u Mhz CPU with %u MB of memory\n", 
				hw.clockMHz, hw.memSize);
		printf("%d torus, %d nodes/pset\n", 
				hw.torus_dimension, hw.sizeOfPset);
		printf("torus size: (%u,%u,%u,%u,%u)\n", 
				hw.Size[0], hw.Size[1],
				hw.Size[2], hw.Size[3],
				hw.Size[4]
		      );
		printf("isTorus?: (%u,%u,%u,%u,%u)\n",
				hw.isTorus[0], hw.isTorus[1], 
				hw.isTorus[2], hw.isTorus[3],
				hw.isTorus[4]
		      );
		iscomplete = 1;
		for(i=0;i<6;i++)
		{
			if(hw.Size[i] > 0 && hw.isTorus[i] != 1)
				iscomplete = 0;
		}
		if(iscomplete)
			printf("Detected a complete torus!\n");
	}
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
}

#endif
