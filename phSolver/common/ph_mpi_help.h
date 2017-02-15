#ifndef PH_MPI_HELP_H
#define PH_MPI_HELP_H

#ifdef __cplusplus
extern "C" {
#endif
int ph_self();
int ph_peers();
int ph_min_int(int val);
int ph_max_int(int val);
long ph_add_long(long val);
double ph_min_double(double val);
double ph_max_double(double val);
double ph_add_double(double val);
#ifdef __cplusplus
}
#endif

#endif
