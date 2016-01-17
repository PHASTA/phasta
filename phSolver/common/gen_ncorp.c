#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include "common_c.h"

#include "FCMangle.h"

//commu_int.f
#define commu_int FortranCInterface_GLOBAL_(commu_int, COMMU_INT)
void commu_int(int* global, int* ilwork, int* n, char* code);

#define gen_ncorp FortranCInterface_GLOBAL_(gen_ncorp, GEN_NCORP)

//  KEJ changed ncorp_t to lcorp_t (used for on proc numbering
//   and introduced gcorp_t for global numbering
typedef int lcorp_t;
#define NCORP_MPI_T MPI_INTEGER
typedef long long int gcorp_t;

static lcorp_t count_owned(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int* num_nodes);
static lcorp_t count_local(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int* num_nodes);


void gen_ncorp(gcorp_t* ncorp, int* ilwork, lcorp_t* nilwork_f, int* num_nodes)
{
	int part;
	int num_parts;
	int i;
	lcorp_t nilwork = *nilwork_f;
	lcorp_t owned;
	lcorp_t local;
	lcorp_t* owner_counts;
	gcorp_t  local_start_id;
	gcorp_t  gid;

	MPI_Comm_rank(MPI_COMM_WORLD, &part);
	MPI_Comm_size(MPI_COMM_WORLD, &num_parts);

	memset(ncorp, 0, sizeof(gcorp_t)*(*num_nodes));
	owned = count_owned(ilwork, nilwork, ncorp, num_nodes);
	local = count_local(ilwork, nilwork, ncorp, num_nodes);
	conpar.iownnodes = owned+local;
#ifdef PRINT_EVERYTHING
	printf("%d: %d local only nodes\n", part, local);
	printf("%d: %d owned nodes\n", part, owned);
#endif
	assert( owned <= *num_nodes );
	assert( owned+local <= *num_nodes );

	owner_counts = (lcorp_t*) malloc(sizeof(lcorp_t)*num_parts);
	memset(owner_counts, 0, sizeof(lcorp_t)*num_parts);
	owner_counts[part] = owned+local;
#ifdef PRINT_EVERYTHING
	for(i=0;i<num_parts;i++)
	{
		printf("%d,", owner_counts[i]);
	}
	printf("\n");
#endif
	MPI_Allgather(MPI_IN_PLACE, 1, NCORP_MPI_T, owner_counts,
		       	1, NCORP_MPI_T, MPI_COMM_WORLD);
#ifdef PRINT_EVERYTHING
	for(i=0;i<num_parts;i++)
	{
		printf("%d,", owner_counts[i]);
	}
	printf("\n");
#endif
	local_start_id=0;
	for(i=0;i<part;i++) //TODO: MPI_Exscan()?
	{
// global so needs long long
		local_start_id += owner_counts[i];
	}
	local_start_id++; //Fortran numbering
#ifdef PRINT_EVERYTHING
	printf("%d: %d\n", part, local_start_id);
#endif
// global so needs long long
	gid = local_start_id;
        if(gid<0) printf("part,gid, %d %ld",part,gid);
        assert(gid>=0);
	for(i=0;i<*num_nodes;i++) //assign owned node's numbers
	{
		//if shared, owned 1
			//if shared, slave -1
			//if local only, 0
		if(ncorp[i] == 1)
		{
// global so needs long long
			ncorp[i]=gid;
                        assert(ncorp[i]>=0);

// global so needs long long
			gid++;
			continue;
		}
		if(ncorp[i] == 0)
		{
			ncorp[i] = gid;
                        assert(ncorp[i]>=0);
			gid++;
			continue;
		}
		if(ncorp[i] == -1)
		{
			ncorp[i] = 0; //commu() adds, so zero slaves
		}

	}
	//char code[] = "out";
	//int ione = 1;
	//commu_int(ncorp, ilwork, &ione, code);

}

static lcorp_t count_local(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int* num_nodes)
{
	int i;
	lcorp_t num_local = 0;
	for(i=0;i<*num_nodes;i++)
	{
		if(ncorp_tmp[i] == 0)
			num_local++; //nodes away from part boundary
		assert(!(ncorp_tmp[i] < -1 || ncorp_tmp[i] > 1));
	}
	return(num_local);
}
static lcorp_t count_owned(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int* num_nodes)
{
	int numtask = ilwork[0];
	int itkbeg = 0; //task offset
	int owned = 0;
	int i,j,k;
	for(i=0;i<numtask;i++)
	{
		int itag = ilwork[itkbeg+1]; //mpi tag
		int iacc = ilwork[itkbeg+2]; //0 for slave, 1 for master
		assert(iacc >= 0 && iacc <= 1);
		int iother = ilwork[itkbeg+3]-1; //other rank (see ctypes.f for off by one)
		int numseg = ilwork[itkbeg+4]; //number of segments
		for(j=0;j<numseg;j++)
		{
			int isgbeg = ilwork[itkbeg+5+(j*2)]; //first idx of seg
			int lenseg = ilwork[itkbeg+6+(j*2)]; //length of seg
			assert(iacc == 0 || iacc == 1);
			if(iacc)
			{
				for(k=0;k<lenseg;k++)
				{
					if(ncorp_tmp[isgbeg-1+k] == 0)
						owned++;
					//make sure we're not both master and slave
					assert(ncorp_tmp[isgbeg-1+k] != -1);
					ncorp_tmp[isgbeg-1+k] = 1;
					assert(isgbeg-1+k < *num_nodes);
				}
				assert(owned <= *num_nodes);
			}
			else
			{
				for(k=0;k<lenseg;k++)
				{
					ncorp_tmp[isgbeg-1+k] = -1;
					assert(isgbeg-1+k < *num_nodes);
				}
			}
			//ncorp_tmp init'd to 0
			//if shared, owned 1
			//if shared, slave -1
			//if local only, 0

			assert(itkbeg+6+(j*2) < nlwork);
		}
		itkbeg+= 4+2*numseg;
	}
	return(owned);
}
