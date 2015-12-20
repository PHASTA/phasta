#ifndef PHSOLVER_PHSTR_H
#define PHSOLVER_PHSTR_H

#ifdef __cplusplus
extern "C" {
#endif
  void phstr_appendInt(char* dest, int v);
  void phstr_appendDbl(char* dest, double src);
  void phstr_appendStr(char* dest, char* src);
#ifdef __cplusplus
}
#endif
#endif

