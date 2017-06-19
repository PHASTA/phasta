#ifndef PHIOSTATS_H
#define PHIOSTATS_H

#include<stdlib.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif
void phastaio_initStats();
void phastaio_addReadBytes(size_t bytes);
void phastaio_addWriteBytes(size_t bytes);
void phastaio_addWriteTime(size_t time);
void phastaio_addReadTime(size_t time);
void phastaio_addOpenTime(size_t time);
void phastaio_addCloseTime(size_t time);
void phastaio_printStats();
#ifdef __cplusplus
}
#endif

#endif
