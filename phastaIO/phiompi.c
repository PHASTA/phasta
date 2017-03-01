#include<mpi.h>
#include<assert.h>
#include"phiompi.h"
size_t phio_ar_sizet(size_t val, int op) {
  size_t res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_UNSIGNED_LONG,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
double phio_ar_dbl(double val, int op) {
  double res = 0;
  int err = MPI_Allreduce(&val,&res,1,MPI_DOUBLE,op,MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);
  return res;
}
size_t phio_min_sizet(size_t val) {
  return phio_ar_sizet(val,MPI_MIN);
}
size_t phio_max_sizet(size_t val) {
  return phio_ar_sizet(val,MPI_MAX);
}
size_t phio_add_sizet(size_t val) {
  return phio_ar_sizet(val,MPI_SUM);
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
void phio_barrier() {
  MPI_Barrier(MPI_COMM_WORLD);
}
