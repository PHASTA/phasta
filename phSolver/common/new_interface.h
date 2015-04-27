#ifndef __NEW_INTERFACE_H__
#define __NEW_INTERFACE_H__

#include <FCMangle.h>

#define igetMinMaxAvg FortranCInterface_GLOBAL_(igetminmaxavg,IGETMINMAXAVG)
#define rgetMinMaxAvg FortranCInterface_GLOBAL_(rgetminmaxavg,RGETMINMAXAVG)
#define print_mesh_stats FortranCInterface_GLOBAL_(print_mesh_stats,PRINT_MESH_STATS)
#define print_mpi_stats FortranCInterface_GLOBAL_(print_mpi_stats,PRINT_MPI_STATS)
#define print_system_stats FortranCInterface_GLOBAL_(print_system_stats,PRINT_SYSTEM_STATS)
#define Write_Restart FortranCInterface_GLOBAL_(write_restart,WRITE_RESTART)
#define Write_Error   FortranCInterface_GLOBAL_(write_error,WRITE_ERROR)
#define Write_Displ   FortranCInterface_GLOBAL_(write_displ,WRITE_DISPL)
#define Write_Field   FortranCInterface_GLOBAL_(write_field,WRITE_FIELD)
//MR CHANGE END
#define Write_PhAvg   FortranCInterface_GLOBAL_(write_phavg,WRITE_PHAVG)
#define Write_PhAvg2  FortranCInterface_GLOBAL_(write_phavg2,WRITE_PHAVG2)
#define Write_d2wall  FortranCInterface_GLOBAL_(write_d2wall,WRITE_D2WALL)
//MR CHANGE END
//NM CHANGE BEGIN
#define Write_Debug   FortranCInterface_GLOBAL_(write_debug,WRITE_DEBUG)
#define POSIX 1
#define SYNCIO 2
#define WRITE_DEBUG_OUTPUT_TYPE POSIX
//NM CHANGE END
#define read_d2wall FortranCInterface_GLOBAL_(read_d2wall,READ_D2WALL)

extern char phasta_iotype[80];
extern int field_flag;
extern int f_descriptor;

void igetMinMaxAvg(int *ivalue, double *stats, int *statRanks);
void rgetMinMaxAvg(double *value, double *stats, int *statRanks);
void print_mesh_stats(void);
void print_mpi_stats(void);
void print_system_stats(double *tcorecp, double *tcorecpscal);

void countfieldstowriterestart();
void
Write_Restart(  int* pid,
                int* stepno,
                int* nshg,
                int* numVars,
                double* array1,
                double* array2 );

void
Write_Error(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 );

void
Write_Displ(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 );
void
Write_Field(  int *pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno);

void
Write_Debug(int *pid,			//ramk of mpi process
			char* filetag,		//Name of the file to write. This will be appended with %d.%d to track step and processor. WARNING: this must be null terminated. 
			char* fieldtag,		//name of field to write data into. WARNING: this must be null terminated. When calling from Fortran, use //char(0) 
			void* array,		//data array to write to file
			char* arraytype,	//arraytype[0] == 'd': double, arraytype[0] == 'i': int
			int* nshg,			//total number of shape functions
			int* numvars,		//number of variables per state (i.e. size of the matrix in dimension 2)
			int* stepno);		//step number

void
Write_Field(  int *pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno);

void
Write_PhAvg2( int* pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              int* iphase,
              int* nphasesincycle,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno);
void
Write_d2wall(   int* pid,
                int* numnp,
                double* array1 );	

void
read_d2wall(  int* pid,
              int* numnp,
              double* array1,
              int* foundd2wall );


#endif //header guard
