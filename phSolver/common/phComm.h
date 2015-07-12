#ifndef PHSOLVER_COMM_H
#define PHSOLVER_COMM_H

#include <FCMangle.h>

#define phcomm_rank FortranCInterface_GLOBAL_(phcomm_rank, PHCOMM_RANK)

#ifdef __cplusplus
extern "C" {
#endif
  int phcomm_rank();
#ifdef __cplusplus
}
#endif

#endif
