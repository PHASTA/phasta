c_serial_test(linkProcsDir-posix
  ln -snf ${CDIR}/2-procs_case-Posix ${CDIR}/2-procs_case)
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
