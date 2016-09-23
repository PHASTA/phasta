set(testLabel "phsolver_compressible_shocktube_${solver}")
add_test(NAME ${testLabel}_posix
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/2-procs_case-Posix
  -DTGTCASEDIR=${CDIR}/2-procs_case
  -DNUMSTART=${CDIR}/numstart.dat
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaC.exe
  -DCOMPARE_EXE=${PHASTA_BINARY_DIR}/bin/checkphasta
  -DIS_SYNCIO=0
  -DNUMPROCS=2
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )

add_test(NAME ${testLabel}_restart-posix
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/2-procs_case-Posix
  -DTGTCASEDIR=${CDIR}/2-procs_case
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaC.exe
  -DCOMPARE_EXE=${PHASTA_BINARY_DIR}/bin/checkphasta
  -DIS_SYNCIO=0
  -DNUMPROCS=2
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )
