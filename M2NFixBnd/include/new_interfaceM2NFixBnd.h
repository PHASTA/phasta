#ifndef __NEW_INTERFACEM2NFIXBND_H__
#define __NEW_INTERFACEM2NFIXBND_H__

#include <FCMangle.h>
#include <mpi.h>

#define Write_M2NFixBnd   FortranCInterface_GLOBAL_(write_m2nfixbnd,WRITE_M2NFIXBND)
#define Write_M2NFixBnd_SolOnly FortranCInterface_GLOBAL_(write_m2nfixbnd_solonly,WRITE_M2NFIXBND_SOLONLY)
#define Write_Restart FortranCInterface_GLOBAL_(write_restart,WRITE_RESTART)
#define Write_Error   FortranCInterface_GLOBAL_(write_error,WRITE_ERROR)
#define Write_Displ   FortranCInterface_GLOBAL_(write_displ,WRITE_DISPL)
#define Write_Field   FortranCInterface_GLOBAL_(write_field,WRITE_FIELD)
#define Write_PhAvg   FortranCInterface_GLOBAL_(write_phavg,WRITE_PHAVG)
#define Write_PhAvg2  FortranCInterface_GLOBAL_(write_phavg2,WRITE_PHAVG2)
#define Write_d2wall  FortranCInterface_GLOBAL_(write_d2wall,WRITE_D2WALL)
#define read_d2wall   FortranCInterface_GLOBAL_(read_d2wall,READ_D2WALL)

extern char phasta_iotype[80];
extern int field_flag;
extern int f_descriptor;

void
Write_M2NFixBnd(int* pid,
                int* stepno,
                int* nshg,
                int* numVars,
                int* ndofybar,
                int* ndoferrors,
                double* array1,
                double* array2,
                double* array3,
                double* array4);

void
Write_M2NFixBnd_SolOnly( int* pid,
                       int* stepno,
                       int* nshg,
                       int* numVars,
                       double* array1 );

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
