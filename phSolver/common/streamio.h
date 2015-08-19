#ifndef STREAMIO_H
#define STREAMIO_H
#include "phstream.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct phio_file* phio_fp;
void streamio_setup_read(phio_fp* f, grstream grs);
void streamio_setup_write(phio_fp* f, rstream rs);
#ifdef __cplusplus
}
#endif
#endif
