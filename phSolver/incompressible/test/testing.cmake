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
