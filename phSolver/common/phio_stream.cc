#include "phio_stream.h"
void stream_openfile(
    const char filename[],
    phio_fp f) {
  stream_fp sf = (stream_fp) f;
  if(sf->mode == 'w' && sf->rs != NULL)
    sf->file = (int*) openRStreamWrite(sf->rs);
  else if(sf->mode == 'r' && sf->grs != NULL)
    sf->file = (int*) openGRStreamRead(sf->grs, filename);
  else
    fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, filename);
}

void stream_closefile(phio_fp f) {
  stream_fp sf = (stream_fp) f;
  fclose((FILE*)sf->file);
}

void stream_constructname(const char* in, char* out) {
  sprintf(out, "%s", in); 
}

