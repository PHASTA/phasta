#ifndef POSIXIO_H
#define POSIXIO_H
#include "phIO.h"
#ifdef __cplusplus
extern "C" {
#endif
void posixio_setup(phio_fp* f, char mode);
#ifdef __cplusplus
}
#endif
#endif
