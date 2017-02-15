#include<mpi.h>
#include<assert.h>
#include"ph_mpi_help.h"
int ph_ar_int(int val, int op) {
  int res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_INT,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
double ph_ar_dbl(double val, int op) {
  double res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_DOUBLE,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
long ph_ar_long(long val, int op) {
  long res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_LONG,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
int ph_min_int(int val) {
  return ph_ar_int(val,MPI_MIN);
}
int ph_max_int(int val) {
  return ph_ar_int(val,MPI_MAX);
}
long ph_add_long(long val) {
  return ph_ar_long(val,MPI_SUM);
}
double ph_min_double(double val) {
  return ph_ar_dbl(val,MPI_MIN);
}
double ph_max_double(double val) {
  return ph_ar_dbl(val,MPI_MAX);
}
double ph_add_double(double val) {
  return ph_ar_dbl(val,MPI_SUM);
}
int ph_self() {
  int self;
  MPI_Comm_rank(MPI_COMM_WORLD, &self);
  return self;
}
int ph_peers() {
  int peers;
  MPI_Comm_size(MPI_COMM_WORLD, &peers);
  return peers;
}
