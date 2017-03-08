#include <phiotimer.h>
void phastaio_time(phastaioTime*) {}
size_t phastaio_time_diff(phastaioTime*, phastaioTime*) {
  return 1;
}
void phastaio_addReadBytes(size_t) {}
void phastaio_addWriteBytes(size_t) {}
void phastaio_addReadTime(size_t) {}
void phastaio_addWriteTime(size_t) {}
void phastaio_setfile(int) {}
void phastaio_addOpenTime(size_t) {}
void phastaio_addCloseTime(size_t) {}
void phastaio_printStats() {}
void phastaio_initStats() {}
