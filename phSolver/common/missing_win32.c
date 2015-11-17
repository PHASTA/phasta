#ifdef intel
#include <stdlib.h>

void SYSTEM(char *cmd) {
    return;
}
void link(char *src, char *dst) {
    return;
}
void  bzero(void* ptr, size_t sz) {
    int i;
    char *cptr;
    cptr = (char*) ptr;
    for (i=0; i < sz; i++) {
        cptr[i]=0;
    }
    return;
}
#else
void dontComplain() {}
#endif
