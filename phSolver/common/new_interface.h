#ifndef __NEW_INTERFACE_H__
#define __NEW_INTERFACE_H__

#include <FCMangle.h>
#include "phIO.h"

#define igetMinMaxAvg FortranCInterface_GLOBAL_(igetminmaxavg,IGETMINMAXAVG)
#define rgetMinMaxAvg FortranCInterface_GLOBAL_(rgetminmaxavg,RGETMINMAXAVG)
#define print_mesh_stats FortranCInterface_GLOBAL_(print_mesh_stats,PRINT_MESH_STATS)
#define print_mpi_stats FortranCInterface_GLOBAL_(print_mpi_stats,PRINT_MPI_STATS)
#define print_system_stats FortranCInterface_GLOBAL_(print_system_stats,PRINT_SYSTEM_STATS)
#define Write_Restart FortranCInterface_GLOBAL_(write_restart,WRITE_RESTART)
#define Write_Error   FortranCInterface_GLOBAL_(write_error,WRITE_ERROR)
#define Write_Field   FortranCInterface_GLOBAL_(write_field,WRITE_FIELD)
#define Write_PhAvg   FortranCInterface_GLOBAL_(write_phavg,WRITE_PHAVG)
#define Write_PhAvg2  FortranCInterface_GLOBAL_(write_phavg2,WRITE_PHAVG2)
#define read_d2wall FortranCInterface_GLOBAL_(read_d2wall,READ_D2WALL)

extern char phasta_iotype[80];
extern int field_flag;
extern phio_fp f_descriptor;

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


#endif
