#ifndef PHSOLVER_COMM_H
#define PHSOLVER_COMM_H

#include <FCMangle.h>

#define phcomm_rank FortranCInterface_GLOBAL_(phcomm_rank, PHCOMM_RANK)
#define phcomm_size FortranCInterface_GLOBAL_(phcomm_size, PHCOMM_SIZE)

#ifdef __cplusplus
extern "C" {
#endif
  int phcomm_rank();
  int phcomm_size();
#ifdef __cplusplus
}
#endif

#endif
