#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include "streamio.h"
#include "phio_stream.h"
#include "phio_posix.h"
#include "phio_base.h"

grstream geomRestartStream;
rstream restartStream;

static struct phio_ops stream_ops = {
  stream_openfile,
  stream_closefile,
  stream_readheader,
  stream_writeheader,
  stream_readdatablock,
  stream_writedatablock,
  stream_constructname
};

void streamio_setup_r(phio_fp* f, RStream* rs, char mode) {
  assert(rs);
  *f = (phio_fp) malloc(sizeof(struct streamio_file));
  stream_fp sf = (stream_fp) *f;
  sf->ops = &stream_ops;
  sf->mode = mode;
  sf->grs = NULL;
  sf->rs = rs;
}

void streamio_setup_read(phio_fp* f, GRStream* grs) {
  *f = (phio_fp) malloc(sizeof(struct streamio_file));
  stream_fp sf = (stream_fp) *f;
  sf->ops = &stream_ops;
  sf->mode = 'r';
  sf->grs = grs;
  sf->rs = NULL;
}

void streamio_setup_write(phio_fp* f, RStream* rs) {
  *f = (phio_fp) malloc(sizeof(struct streamio_file));
  stream_fp sf = (stream_fp) *f;
  sf->ops = &stream_ops;
  sf->mode = 'w';
  sf->grs = NULL;
  sf->rs = rs;
}

void streamio_set_gr(grstream grs) {
  geomRestartStream = grs;
}

grstream streamio_get_gr() {
  return geomRestartStream;
}

void streamio_set_r(rstream rs) {
  restartStream = rs;
}

rstream streamio_get_r() {
  return restartStream;
}
