#include "phstream.h"
void* fail(const char* f) {
  fprintf(stderr,
    "ERROR: function %s is disabled - compile with chefPhasta\n", f);
  return NULL;
}
rstream makeRStream() {
  return (rstream)fail(__func__);
}
void clearRStream(rstream rs) {
  fail(__func__);
}
void destroyRStream(rstream rs) {
  fail(__func__);
}

grstream makeGRStream() {
  return (grstream)fail(__func__);
}
void clearGRStream(grstream grs) {
  fail(__func__);
}
void destroyGRStream(grstream grs) {
  fail(__func__);
}

FILE* openRStreamRead(rstream rs) {
  return (FILE*)fail(__func__);
}
FILE* openRStreamWrite(rstream rs) {
  return (FILE*)fail(__func__);
}

FILE* openGRStreamRead(grstream grs, const char* named) {
  return (FILE*)fail(__func__);
}
FILE* openGRStreamWrite(grstream grs, const char* named) {
  return (FILE*)fail(__func__);
}

void attachRStream(grstream grs, rstream rs) {
  fail(__func__);
}
