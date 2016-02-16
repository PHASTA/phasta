c_serial_test(linkProcsDir-sync
  ln -snf ${CDIR}/2-procs_case-SyncIO-1 ${CDIR}/2-procs_case)
if(HAS_VALGRIND)
  c_serial_test(resetNumStartValgrind-sync
    cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
  c_parallel_test(valgrind-sync 2 ${CDIR}
    valgrind --leak-check=yes --log-file=cSyncValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
endif(HAS_VALGRIND)
c_serial_test(resetNumStart-sync
  cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
c_parallel_test( sync 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
set(compareArgs
  ${CDIR}/2-procs_case-SyncIO-1/
  ${CDIR}/2-procs_case-SyncIO-1_ref/
  1 1e-6)
c_parallel_test(compare-sync 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
c_parallel_test(restart-sync 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
c_parallel_test(compareRestart-sync 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})

c_serial_test(linkProcsDir-posix
  ln -snf ${CDIR}/2-procs_case-Posix ${CDIR}/2-procs_case)
if(HAS_VALGRIND)
  c_serial_test(resetNumStartValgrind-posix
    cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
  c_parallel_test(compressibleValgrind-posix 2 ${CDIR}
    valgrind --leak-check=yes --log-file=cPosixValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
endif(HAS_VALGRIND)
c_serial_test(resetNumStart-posix
  cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
c_parallel_test(posix 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
set(compareArgs
  ${CDIR}/2-procs_case-Posix/
  ${CDIR}/2-procs_case-Posix_ref/
  0 1e-6)
c_parallel_test(compare-posix 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
c_parallel_test(restart-posix 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
c_parallel_test(compareRestart-posix 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
c_serial_test(unlinkProcsDir-compressible rm ${CDIR}/2-procs_case)
