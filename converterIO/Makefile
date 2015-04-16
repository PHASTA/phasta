cxx := mpicxx
cc := mpicc

#mpi_path := /usr/local/mpich2-1.2/
#mpi_libs := -lmpich


all: converterO2N converterN2O

converterO2N: converterO2N.o  phastaIO.o
	$(cxx)  converterO2N.o phastaIO.o -g -o converterO2N

converterN2O: converterN2O.o  phastaIO.o
	$(cxx)  converterN2O.o phastaIO.o -g -o converterN2O

phastaIO.o:
	$(cxx)  -g -DMPICH_IGNORE_CXX_SEEK -c phastaIO.cc -o phastaIO.o

converterO2N.o:
	$(cxx)  -g -DMPICH_IGNORE_CXX_SEEK -c converterO2N.cc -o converterO2N.o


converterN2O.o:
	 $(cxx)  -g -DMPICH_IGNORE_CXX_SEEK -c converterN2O.cc -o converterN2O.o

clean:
	rm -rf *.o converterO2N converterN2O
