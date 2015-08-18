#ifndef STREAMIO_H
#define STREAMIO_H
#ifdef __cplusplus
extern "C" {
#endif
typedef struct phio_file* phio_fp;
typedef struct Stream stream;
void streamio_setup(stream* s, phio_fp* f, char mode);
#ifdef __cplusplus
}
#endif
#endif
