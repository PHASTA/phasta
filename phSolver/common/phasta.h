#ifndef PHASTA_H_
#define PHASTA_H_
struct RStream;
struct GRStream;
int phasta(int argc, char**argv);
int phasta(int argc, char**argv, GRStream* in);
#endif
