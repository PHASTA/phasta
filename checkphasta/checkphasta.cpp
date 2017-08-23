#ifndef __xlC__
#define _GNU_SOURCE
#endif
#define _WITH_DPRINTF
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

#include "syncio.h"
#include "posixio.h"
#include "phIO.h"

char read_solution(double** solutiono, int* size, int* nshgo, int* ndofo,
		int nump, int rank, int timestep, int nSyncFiles, char* casedir);
std::set<int>* find_timesteps(char* casedir, int nSyncFiles);
double compare_solution(char* lpath, char* rpath, 
    int timestep, int nump, int nSyncFiles);
char* getRestartName(int nSyncFiles);

int main(int argc, char** argv)
{
       	int rank;
	int size;
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
        
	if(argc != 5) {
          fprintf(stderr, "argc %d\n", argc);
          fprintf(stderr, 
              "Usage: %s <left> <right> <numSyncFiles> <tolerance>\n"
              "where <left> and <right> are different"
              "N-procs_case directories\n", argv[0]);
          MPI_Finalize();
          return 1;
        }
	char* lpath = argv[1];
	char* rpath = argv[2];
        int nSyncFiles = atoi(argv[3]);
        double tolerance = atof(argv[4]);

	int ndof;
	int nshg;
	int solsize;
	double* solution;

	std::set<int>* l_timesteps = find_timesteps(lpath, nSyncFiles);
	std::set<int>* r_timesteps = find_timesteps(rpath, nSyncFiles);
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
	read_solution(&solution, &solsize, &nshg, &ndof, 
            size, rank, 0, numSyncFiles, "./");
	printf("nshg: %d, ndof: %d\n", nshg, ndof);
	assert(solsize == ndof*nshg);
#endif
	double maxerror = 0.0;
	double error;
	double gblmaxerror;
	for(std::set<int>::iterator i = timesteps_to_check->begin();
			i!=timesteps_to_check->end();i++)
	{
		error = compare_solution(lpath, rpath, *i, size, nSyncFiles);
		if(error>maxerror) maxerror = error;
	}
        delete timesteps_to_check;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&maxerror, &gblmaxerror, 1, MPI_DOUBLE, MPI_MAX, 0,
		MPI_COMM_WORLD);
	if(rank == 0) printf("Maximum difference across all timesteps: %e\n", 
			gblmaxerror);
	MPI_Finalize();
        return (gblmaxerror > tolerance);
}
double compare_solution(char* lpath, char* rpath, int timestep, int nump, int nSyncFiles)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double* lsol;
	double* rsol;
	int lsize;
	int rsize;

	read_solution(&lsol, &lsize, NULL, NULL, 
            nump, rank, timestep, nSyncFiles, lpath);
	read_solution(&rsol, &rsize, NULL, NULL, 
            nump, rank, timestep, nSyncFiles, rpath);

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
		int nump, int rank, int timestep, int nSyncFiles, char* casedir)
{
	int iarray[10];
        const char* iformat = "binary";
        int ithree=3;
        int igeombc;
        char* fn;
	int nshg;
	int ndof;
	double* solution;
        phio_fp fp;
        if( nSyncFiles == 0 )
          posixio_setup(&fp, 'r');
        else if( nSyncFiles > 0 )
          syncio_setup_read(nSyncFiles, &fp);
        char rname[1024];
        phio_constructName(fp,"restart",rname);
        asprintf(&fn,"%s/%s%d.",casedir,rname,timestep);
        phio_openfile(fn, fp);

	phio_readheader(fp, "solution", (void*) iarray, &ithree, "integer", iformat);
	nshg = iarray[0];
	ndof = iarray[1];
	if(size != NULL)
		*size = nshg*ndof;
	solution = (double*) malloc(sizeof(double)*nshg*ndof);
	phio_readdatablock(fp, "solution", solution, size, "double", iformat);
	phio_closefile(fp);
	if(solutiono != NULL)
		*solutiono = solution;
	if(nshgo != NULL)
		*nshgo = nshg;
	if(ndofo != NULL)
		*ndofo = ndof;
	free(fn);
	return(0);
}

std::set<int>* find_timesteps(char* casedir, int nSyncFiles)
{
	char* path;
	DIR* dir;
	struct dirent* d;
	int part, ts;
	std::set<int>* step_list = new std::set<int>;

        asprintf(&path, "%s", casedir);
	dir = opendir(path);
	if(!dir)
	{
		perror("Error opening case: "); 
		MPI_Abort(MPI_COMM_WORLD,1);
	}
        char* rname = getRestartName(nSyncFiles);
        char* fmt;
        asprintf(&fmt, "%s.%%d.%%d", rname);
	while((d=readdir(dir)))
	{
		if(sscanf(d->d_name, fmt, &ts, &part)==2)
		{
			step_list->insert(ts);
		}
	}
        free(rname);
        free(fmt);
	free(path);
        closedir(dir);
	return(step_list);
}

char* getRestartName(int nSyncFiles) {
  char* f;
  if(0 == nSyncFiles)
    asprintf(&f, "restart");
  else if(nSyncFiles > 0)
    asprintf(&f, "restart-dat");
  else {
    fprintf(stderr, 
        "ERROR: the number of sync-io files must be"
        "greater than or equal to zero\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    return NULL;
  }
  return f;
}
