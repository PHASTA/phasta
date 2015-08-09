#ifndef POSIXIO_H
#define POSIXIO_H
#ifdef __cplusplus
extern "C" {
#endif
typedef struct phio_file* phio_fp;
void posixio_setup(phio_fp* f, char mode);
#ifdef __cplusplus
}
#endif
#endif
