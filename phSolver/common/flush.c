#include <FCMangle.h>
#define flush FortranCInterface_GLOBAL(flush, FLUSH)

void flush(int* junk ){ return; }

