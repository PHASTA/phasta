#ifndef __PHIOTMRC_H__
#define __PHIOTMRC_H__

#include<stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif
double phiotmrc (void);

#ifdef __INTEL_COMPILER
typedef size_t phioTime;
#else
typedef struct timespec phioTime;
#endif
void phastaio_time(phioTime* t);
size_t phastaio_time_diff(phioTime* start, phioTime* end);
void phastaio_setCyclesPerMicroSec();
#ifdef __cplusplus
}
#endif

#endif // __PHIOTMRC_H__
