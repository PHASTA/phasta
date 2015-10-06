#ifndef SYNCIO_H
#define SYNCIO_H
#include "phIO.h"
#ifdef __cplusplus
extern "C" {
#endif
void syncio_setup_read(int nfiles, phio_fp* f);
void syncio_setup_write(int nfiles, int nfields, int nppf, phio_fp* f);
#ifdef __cplusplus
}
#endif
#endif
