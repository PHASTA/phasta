#include <sys/stat.h>
#include <sys/types.h>

int ph_mkdir(const char* path) {
  int err = mkdir(path, S_IRWXU);
  return err;
}
