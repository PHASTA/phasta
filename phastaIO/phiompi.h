#ifndef PHIOMPI_H
#define PHIOMPI_H

#include<stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif
int phio_self();
int phio_peers();
void phio_barrier();
size_t phio_min_sizet(size_t val);
size_t phio_max_sizet(size_t val);
size_t phio_add_sizet(size_t val);
double phio_min_double(double val);
double phio_max_double(double val);
double phio_add_double(double val);
#ifdef __cplusplus
}
#endif

#endif
