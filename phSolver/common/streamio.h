#ifndef STREAMIO_H
#define STREAMIO_H
#include "phstream.h"
#include "phIO.h"
#ifdef __cplusplus
extern "C" {
#endif
void streamio_setup_r(phio_fp* f, rstream rs, char mode);
void streamio_setup_read(phio_fp* f, grstream grs);
void streamio_setup_write(phio_fp* f, rstream rs);
void streamio_set_gr(grstream grs);
grstream streamio_get_gr();
void streamio_set_r(rstream rs);
rstream streamio_get_r();
#ifdef __cplusplus
}
#endif
#endif
