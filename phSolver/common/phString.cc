#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <sstream>

#ifdef __cplusplus
extern "C" {
#endif

void phstr_appendInt(char* dest, int v) {
  std::stringstream ss;
  ss << dest << v;
  std::string s = ss.str();
  strcpy(dest, s.c_str());
}

void phstr_appendDbl(char* dest, double v) {
  std::stringstream ss;
  ss << dest << v;
  std::string s = ss.str();
  strcpy(dest, s.c_str());
}

void phstr_appendStr(char* dest, char* src) {
  std::stringstream ss;
  ss << dest << src;
  std::string s = ss.str();
  strcpy(dest, s.c_str());
}

#ifdef __cplusplus
}
#endif
