#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> 
#include <set>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include "phastaIO.h"

//Provided by phastaIO
void Gather_Headers( int* fileDescriptor, std::vector< std::string >& headers );

char read_solution(double** solutiono, int* size, int* nshgo, int* ndofo,
		int nump, int rank, int timestep, char* casedir);

std::set<int>* find_timesteps(char* casedir, int nump);
double compare_solution(char* lpath, char* rpath, int timestep, int nump);

int main(int argc, char** argv)
{
	int rank;
	int size;
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	assert(argc>2);
	char* lpath = argv[1];
	char* rpath = argv[2];

	int ndof;
	int nshg;
	int solsize;
	double* solution;

	std::set<int>* l_timesteps = find_timesteps(lpath, size);
	std::set<int>* r_timesteps = find_timesteps(rpath, size);
	std::set<int>* timesteps_to_check = new std::set<int>;
	std::set_intersection(l_timesteps->begin(), l_timesteps->end(),
			r_timesteps->begin(), r_timesteps->end(),
			std::inserter(*timesteps_to_check, timesteps_to_check->begin()));
        delete l_timesteps;
        delete r_timesteps;
	if(rank == 0)
		printf("Found %d common timesteps\n",
			       	timesteps_to_check->size());
#ifdef DBGONLY
	read_solution(&solution, &solsize, &nshg, &ndof, size, rank, 0, "./");
	printf("nshg: %d, ndof: %d\n", nshg, ndof);
	assert(solsize == ndof*nshg);
#endif
	double maxerror = 0.0;
	double error;
	double gblmaxerror;
	for(std::set<int>::iterator i = timesteps_to_check->begin();
			i!=timesteps_to_check->end();i++)
	{
		error = compare_solution(lpath, rpath, *i, size);
		if(error>maxerror) maxerror = error;
	}
        delete timesteps_to_check;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&maxerror, &gblmaxerror, 1, MPI_DOUBLE, MPI_MAX, 0,
		MPI_COMM_WORLD);
	if(rank == 0) printf("Maximum difference across all timesteps: %e\n", 
			gblmaxerror);
	MPI_Finalize();
}
double compare_solution(char* lpath, char* rpath, int timestep, int nump)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double* lsol;
	double* rsol;
	int lsize;
	int rsize;

	read_solution(&lsol, &lsize, NULL, NULL, nump, rank, timestep, lpath);
	read_solution(&rsol, &rsize, NULL, NULL, nump, rank, timestep, rpath);

	double maxdiff=0.0;
	double gblmaxdiff;
	if(lsize != rsize)
	{
		printf("Error: Solution sizes different: %d, %d\n", 
				lsize, rsize);
		assert(lsize == rsize);
	}
	for(int i=0;i<lsize;i++)
	{
		double diff = fabs(lsol[i]-rsol[i]);
		if(diff > maxdiff) maxdiff = diff;
	}
        free(lsol);
        free(rsol);
	MPI_Reduce(&maxdiff, &gblmaxdiff, 1, MPI_DOUBLE, MPI_MAX, 0,
		       	MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD); //TODO: debugging only
	if(rank == 0)
	printf("Timestep: %d, maximum difference: %e\n", timestep, gblmaxdiff);
	return gblmaxdiff;

}
char read_solution(double** solutiono, int* size, int* nshgo, int* ndofo,
		int nump, int rank, int timestep, char* casedir)
{
	int iarray[10];
        const char* iformat = "binary";
        int ithree=3;
        int igeombc;
        char* fn;
	int nshg;
	int ndof;
	double* solution;
	if(nump == 1)
		asprintf(&fn, "%s/restart.%d.%d",
				casedir,timestep,rank+1);
	else
		asprintf(&fn, "%s/%d-procs_case/restart.%d.%d",
				casedir, nump,timestep,rank+1);
        openfile(fn, "read", &igeombc);
	//TODO: error handle
	readheader(&igeombc, "solution", (void*) iarray, &ithree, "integer", iformat);
	nshg = iarray[0];
	ndof = iarray[1];
	if(size != NULL)
		*size = nshg*ndof;
	solution = (double*) malloc(sizeof(double)*nshg*ndof);
	readdatablock(&igeombc, "solution", solution, size, "double", iformat);
	closefile(&igeombc, "read");
	if(solutiono != NULL)
		*solutiono = solution;
	if(nshgo != NULL)
		*nshgo = nshg;
	if(ndofo != NULL)
		*ndofo = ndof;
	free(fn);
	return(0);
}

std::set<int>* find_timesteps(char* casedir, int nump)
{
	char* path;
	char* fullpath;
	DIR* dir;
	struct dirent* d;
	int part, ts;
	std::set<int>* step_list = new std::set<int>;

	if(nump == 1)
		asprintf(&path, "%s", casedir);
	else
		asprintf(&path, "%s/%d-procs_case", casedir, nump);
	dir = opendir(path);
	if(!dir)
	{
		perror("Error opening case: "); 
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	while((d=readdir(dir)))
	{
		asprintf(&fullpath, "%s/%s", path, d->d_name);
		if(sscanf(d->d_name, "restart.%d.%d", &ts, &part)==2)
		{
			step_list->insert(ts);
		}
	}
	return(step_list);
	free(path);
}

void check_ilwork(int* ilwork, int nlwork)
{
	int numtask = ilwork[0];
	int itkbeg = 0; //task offset
	printf("%d tasks\n", numtask);

	for(int i=0;i<numtask;i++)
	{
		int itag = ilwork[itkbeg+1]; //mpi tag
		int iacc = ilwork[itkbeg+2]; //0 for slave, 1 for master
		int iother = ilwork[itkbeg+3]-1; //other rank (see ctypes.f for off by one)
		int numseg = ilwork[itkbeg+4]; //number of segments
		printf("Comm with rank: %d\n", iother);
		for(int j=0;j<numseg;j++)
		{
			int isgbeg = ilwork[itkbeg+5+(j*2)]; //first idx of seg
			int lenseg = ilwork[itkbeg+6+(j*2)]; //length of seg

			printf("isgbeg: %d, len: %d\n", isgbeg, lenseg);
			assert(itkbeg+6+(j*2) < nlwork);
		}
		itkbeg+= 4+2*numseg;
	}
}

void read_ilwork()
{
	int rank;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int iarray[10];
	const char* iformat = "binary";
	int ione=1;
	int itwo=2;
	int igeombc;
	int nilwork=0;
	int tmp;
	int* ilwork;
	char* fn;
	asprintf(&fn, "%d-procs_case/geombc.dat.%d",size,rank+1);
	openfile(fn, "read", &igeombc);

	std::vector<std::string> headers;
	Gather_Headers(&igeombc, headers);

	readheader(&igeombc, "size of ilwork array", (void*) &nilwork, &ione, 
			"integer", iformat);
	assert(nilwork > 0);
	readheader(&igeombc, "ilwork", (void*) &tmp, &ione, 
			"integer", iformat);
	assert(tmp == nilwork);

	ilwork = (int*) malloc(sizeof(int)*nilwork);

	readdatablock(&igeombc, "ilwork", ilwork, &nilwork, "integer", iformat);
	closefile(&igeombc, "read");
}
