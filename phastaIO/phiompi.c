#include<mpi.h>
#include<assert.h>
#include"phiompi.h"
int phio_ar_int(int val, int op) {
  int res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_INT,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
double phio_ar_dbl(double val, int op) {
  double res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_DOUBLE,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
long phio_ar_long(long val, int op) {
  long res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_LONG,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
int phio_min_int(int val) {
  return phio_ar_int(val,MPI_MIN);
}
int phio_max_int(int val) {
  return phio_ar_int(val,MPI_MAX);
}
long phio_add_long(long val) {
  return phio_ar_long(val,MPI_SUM);
}
double phio_min_double(double val) {
  return phio_ar_dbl(val,MPI_MIN);
}
double phio_max_double(double val) {
  return phio_ar_dbl(val,MPI_MAX);
}
double phio_add_double(double val) {
  return phio_ar_dbl(val,MPI_SUM);
}
int phio_self() {
  int self;
  MPI_Comm_rank(MPI_COMM_WORLD, &self);
  return self;
}
int phio_peers() {
  int peers;
  MPI_Comm_size(MPI_COMM_WORLD, &peers);
  return peers;
}
