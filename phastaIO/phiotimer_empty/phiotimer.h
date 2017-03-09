#ifndef PHIOTIMER_EMPTY_H
#define PHIOTIMER_EMPTY_H

#include<stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PHASTAIO_OPENTIME(cmd) cmd
#define PHASTAIO_CLOSETIME(cmd) cmd
#define PHASTAIO_READTIME(cmd,ignored) cmd
#define PHASTAIO_WRITETIME(cmd,ignored) cmd

enum phastaio_file { GEOMBC_READ, RESTART_READ, RESTART_WRITE };

typedef int phastaioTime;
struct phastaio_stats;
void phastaio_time(phastaioTime*);
size_t phastaio_time_diff(phastaioTime*, phastaioTime*);
void phastaio_addReadBytes(size_t);
void phastaio_addWriteBytes(size_t);
void phastaio_addReadTime(size_t);
void phastaio_addWriteTime(size_t);
void phastaio_setfile(int);
void phastaio_addOpenTime(size_t);
void phastaio_addCloseTime(size_t);
void phastaio_printStats();
void phastaio_initStats();

#ifdef __cplusplus
}
#endif

#endif
