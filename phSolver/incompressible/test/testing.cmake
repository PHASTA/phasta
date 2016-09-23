set(testLabel "phsolver_incompressible")
add_test(NAME ${testLabel}_sync
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-SyncIO-2
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DNUMSTART=${CDIR}/numstart.dat
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )
if(HAS_VALGRIND)
  ic_serial_test(resetNumStartValgrind-sync
    cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
  ic_parallel_test(valgrind-sync 4 ${CDIR}
    valgrind --leak-check=yes --log-file=icSyncValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
endif(HAS_VALGRIND)
set(compareArgs
  ${CDIR}/4-procs_case-SyncIO-2/
  ${CDIR}/4-procs_case-SyncIO-2_ref/
  2 1e-6)
ic_parallel_test(compare-sync 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})

add_test(NAME ${testLabel}_restart-sync
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-SyncIO-2
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )
ic_parallel_test(compareRestart-sync 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})

add_test(NAME ${testLabel}_posix
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-Posix
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DNUMSTART=${CDIR}/numstart.dat
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )

if(HAS_VALGRIND)
  ic_serial_test(resetNumStartValgrind-posix
    cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
  ic_parallel_test(valgrind-posix 4 ${CDIR}
    valgrind --leak-check=yes --log-file=icPosixValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
endif(HAS_VALGRIND)

set(compareArgs
  ${CDIR}/4-procs_case-Posix/
  ${CDIR}/4-procs_case-Posix_ref/
  0 1e-6)
ic_parallel_test(compare-posix 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})

add_test(NAME ${testLabel}_restart-posix
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-Posix
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )
ic_parallel_test(compareRestart-posix 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
