#include <stdio.h>
#include <stdlib.h>
#include "streamio.h"
#include "phio_stream.h"
#include "phio_posix.h"
#include "phio_base.h"

static struct phio_ops stream_ops = {
  stream_openfile_write,
  stream_closefile,
  posix_readheader,
  posix_writeheader,
  posix_readdatablock,
  posix_writedatablock,
  stream_constructname
};

void streamio_setup(phio_fp*, GRStream grs) {
  *f = (phio_fp) malloc(sizeof(struct streamio_file));
  stream_fp sf = (stream_fp) *f;
  f->ops = &stream_ops;
  f->file = (int*) malloc(sizeof(int*));
  f->mode = mode;
  sf->grs = grs;
}

void streamio_setup(phio_fp*, RStream rs) {
  *f = (phio_fp) malloc(sizeof(struct streamio_file));
  stream_fp sf = (stream_fp) *f;
  f->ops = &stream_ops;
  f->file = (int*) malloc(sizeof(int*));
  f->mode = mode;
  sf->rs = rs;
}
