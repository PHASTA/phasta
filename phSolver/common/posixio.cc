#include <stdio.h>
#include <stdlib.h>
#include "posixio.h"
#include "phio_posix.h"
#include "phio_base.h"

static struct phio_ops posix_ops = {
  posix_openfile,
  posix_closefile,
  posix_readheader,
  posix_writeheader,
  posix_readdatablock,
  posix_writedatablock,
  posix_constructname
};

void posixio_setup(phio_fp* f, char mode) {
  *f = (phio_fp) malloc(sizeof(struct phio_file));
  (*f)->ops = &posix_ops;
  if(mode != 'r' && mode != 'w') {
    fprintf(stderr, "ERROR unsupported file mode in %s on line %d"
        "... exiting", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  (*f)->file = (int*) malloc(sizeof(int*));
  (*f)->mode = mode;
}
