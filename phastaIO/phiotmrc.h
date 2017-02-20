#ifndef __PHIOTMRC_H__
#define __PHIOTMRC_H__

#include<stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif
double phiotmrc (void);

typedef size_t phioTime;
void phastaio_time(phioTime* t);
size_t phastaio_time_diff(phioTime* start, phioTime* end);
void phastaio_setCyclesPerMicroSec();
#ifdef __cplusplus
}
#endif

#endif // __PHIOTMRC_H__
