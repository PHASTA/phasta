#include <stdio.h>
#include <stdlib.h>
#include "syncio.h"
#include "phio_sync.h"
#include "phio_base.h"

static struct phio_ops sync_ops_write = {
  sync_openfile_write,
  sync_closefile,
  sync_readheader,
  sync_writeheader,
  sync_readdatablock,
  sync_writedatablock,
  sync_constructname
};

static struct phio_ops sync_ops_read = {
  sync_openfile_read,
  sync_closefile,
  sync_readheader,
  sync_writeheader,
  sync_readdatablock,
  sync_writedatablock,
  sync_constructname
};

void init(sync_fp f, char mode) {
  if(mode == 'w')
    f->ops = &sync_ops_write;
  else if(mode == 'r')
    f->ops = &sync_ops_read;
  else {
    fprintf(stderr, "ERROR unsupported file mode in %s on line %d"
        "... exiting", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  f->file = (int*) malloc(sizeof(int*));
  f->mode = mode;
}

void syncio_setup_read(int nfiles, phio_fp* f) {
  *f = (phio_fp) malloc(sizeof(struct syncio_file));
  sync_fp sf = (sync_fp) *f;
  init(sf, 'r');
  sf->nfiles = nfiles;
  sf->nfields = 0;
  sf->nppf = 0;
}

void syncio_setup_write(int nfiles, int nfields, int nppf, phio_fp* f) {
  *f = (phio_fp) malloc(sizeof(struct syncio_file));
  sync_fp sf = (sync_fp) *f;
  init(sf, 'w');
  sf->nfiles = nfiles;
  sf->nfields = nfields;
  sf->nppf = nppf;
}
