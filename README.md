#build and test
  wget www.scorec.rpi.edu/~cwsmith/phastaCases.tar.gz .
  tar xzf phastaCases.tar.gz

  cmake \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Debug \
  -DPHASTA_INCOMPRESSIBLE=ON \
  -DPHASTA_COMPRESSIBLE=ON \
  -DACUSOLVE_LIB=/path/to/libles.a \
  -DCASES=/path/to/phastaCases/ \
  ..

  make

  ctest
