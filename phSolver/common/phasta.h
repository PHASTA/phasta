#ifndef PHASTA_H_
#define PHASTA_H_
#include "Input.h"
struct RStream;
struct GRStream;
int phasta(int argc, char**argv);
int phasta(phSolver::Input& ctrl);
int phasta(phSolver::Input& ctrl, GRStream* in);
int phasta(phSolver::Input& ctrl, RStream* out);
int phasta(phSolver::Input& ctrl, GRStream* in, RStream* out);
#endif
