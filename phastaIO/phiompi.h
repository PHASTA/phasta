#ifndef PHIOMPI_H
#define PHIOMPI_H

#ifdef __cplusplus
extern "C" {
#endif
int phio_self();
int phio_peers();
int phio_min_int(int val);
int phio_max_int(int val);
long phio_add_long(long val);
double phio_min_double(double val);
double phio_max_double(double val);
double phio_add_double(double val);
#ifdef __cplusplus
}
#endif

#endif
