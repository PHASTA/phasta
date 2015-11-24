#ifndef COMMON_TMRC_H
#define COMMON_TMRC_H

#include <FCMangle.h>
#define TMRC FortranCInterface_GLOBAL_(tmrc, TMRC)

#ifdef __cplusplus
extern "C" {
#endif

double TMRC(void);

#ifdef __cplusplus
}
#endif

#endif
