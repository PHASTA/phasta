/*  Primary interface for the Phasta Binary read and write routines these*/
/*  functions are 'C' callable.( All arguments have been kept as pointers to*/
/*  facilitate calling from Fortran )*/
/*  Michel Rasquin  Spring 2012, inspired from phastaIO.h*/
#ifndef _SETSYNCIOPARAM_H_
#define _SETSYNCIOPARAM_H_

#include <FCMangle.h>

#define detectd2wallfiles FortranCInterface_GLOBAL_(detectd2wallfiles,DETECTD2WALLFILES)


#if defined (__cplusplus)
extern "C" {
#endif

  void detectd2wallfiles(int* numd2wallfiles);

#ifdef __cplusplus
} // end of extern "C".

#endif // __cplusplus

#endif // _SETSYNCIOPARAM_H_
