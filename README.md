#build and test

    wget https://fluid.colorado.edu/~kjansen/PHASTA/phastaChefTests.tar.gz
    tar xzf phastaChefTests.tar.gz # use for CASES path below
    
Note, the following builds only the native compressible solver.  There are options that can be turned on to utilize PETSc for compressible.  Likewise there are options to build the incompressible solver with SVLS. Finally CMAKE can configure alternate compilers and optimization choices. 

    cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_BUILD_TYPE=Debug \
    -DPHASTA_INCOMPRESSIBLE=OFF \
    -DPHASTA_COMPRESSIBLE=ON \
    -DPHASTA_USE_SVLS=OFF \
    -DPHASTA_USE_PETSC=OFF \    
    -DPHASTA_TESTING=ON \
    -DCASES=/path/to/phastaCases/ \
    ..

    make

    ctest
