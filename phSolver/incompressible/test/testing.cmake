macro(ic_parallel_test name procs dir exe)
  set(tname incompressible_${name})
  add_test(
    NAME ${tname}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${procs} ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_incompressible")
endmacro(ic_parallel_test)

macro(ic_serial_test name exe)
  set(tname incompressible_${name})
  add_test( NAME ${tname} COMMAND ${exe} ${ARGN} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_incompressible")
endmacro(ic_serial_test)

set(CDIR ${CASES}/incompressible)
ic_serial_test(copyInpCfg
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
ic_serial_test(linkProcsDir-sync
  ln -snf ${CDIR}/4-procs_case-SyncIO-2 ${CDIR}/4-procs_case)
if(HAS_VALGRIND)
  ic_serial_test(resetNumStartValgrind-sync
    cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
  ic_parallel_test(valgrind-sync 4 ${CDIR}
    valgrind --leak-check=yes --log-file=icSyncValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
endif(HAS_VALGRIND)
ic_serial_test(resetNumStart-sync
  cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
ic_parallel_test(sync 4 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
set(compareArgs
  ${CDIR}/4-procs_case-SyncIO-2/ 
  ${CDIR}/4-procs_case-SyncIO-2_ref/ 
  2 1e-6)
ic_parallel_test(compare-sync 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
ic_parallel_test(restart-sync 4 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
ic_parallel_test(compareRestart-sync 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})

ic_serial_test(linkProcsDir-posix
  ln -snf ${CDIR}/4-procs_case-Posix ${CDIR}/4-procs_case)
if(HAS_VALGRIND)
  ic_serial_test(resetNumStartValgrind-posix
    cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
  ic_parallel_test(valgrind-posix 4 ${CDIR}
    valgrind --leak-check=yes --log-file=icPosixValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
endif(HAS_VALGRIND)
ic_serial_test(resetNumStart-posix
  cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
ic_parallel_test(posix 4 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
set(compareArgs 
  ${CDIR}/4-procs_case-Posix/ 
  ${CDIR}/4-procs_case-Posix_ref/ 
  0 1e-6)
ic_parallel_test(compare-posix 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
ic_parallel_test(restart-posix 4 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaIC.exe)
ic_parallel_test(compareRestart-posix 4 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
ic_serial_test(unlinkProcsDir rm ${CDIR}/4-procs_case)
