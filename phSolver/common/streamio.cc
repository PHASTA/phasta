#include <stdio.h>
#include <stdlib.h>
#include "streamio.h"
#include "phio_stream.h"
#include "phio_posix.h"
#include "phio_base.h"

static struct phio_ops stream_ops = {
  stream_openfile,
  stream_closefile,
  posix_readheader,
  posix_writeheader,
  posix_readdatablock,
  posix_writedatablock,
  stream_constructname
};

void streamio_setup_read(phio_fp* f, GRStream* grs) {
  *f = (phio_fp) malloc(sizeof(struct streamio_file));
  stream_fp sf = (stream_fp) *f;
  sf->ops = &stream_ops;
  sf->file = (int*) malloc(sizeof(int*));
  sf->mode = 'r';
  sf->grs = grs;
  sf->rs = NULL;
}

void streamio_setup_write(phio_fp* f, RStream* rs) {
  *f = (phio_fp) malloc(sizeof(struct streamio_file));
  stream_fp sf = (stream_fp) *f;
  sf->ops = &stream_ops;
  sf->file = (int*) malloc(sizeof(int*));
  sf->mode = 'w';
  sf->grs = NULL;
  sf->rs = rs;
}
