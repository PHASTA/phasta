set(CDIR ${CASES}/compressible)
add_test(copyInpCfg 
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

add_test(linkProcsDir-sync
  ln -snf ${CDIR}/2-procs_case-SyncIO-1 ${CDIR}/2-procs_case)
add_test(compressibleResetNumStart-sync
  cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
add_test(
  NAME compressible-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${PHASTA_BINARY_DIR}/bin/phastaC.exe
  WORKING_DIRECTORY ${CDIR}
)
set(cmd 
  ${PHASTA_BINARY_DIR}/bin/checkphasta 
  ${CDIR}/2-procs_case-SyncIO-1/ 
  ${CDIR}/2-procs_case-SyncIO-1_ref/ 
  1 1e-6)
add_test(
  NAME compareCompressible-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME compressibleRestart-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${PHASTA_BINARY_DIR}/bin/phastaC.exe
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME compareCompressibleRestart-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)

add_test(linkProcsDir-posix
  ln -snf ${CDIR}/2-procs_case-Posix ${CDIR}/2-procs_case)
add_test(compressibleResetNumStart-posix
  cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
add_test(
  NAME compressible-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${PHASTA_BINARY_DIR}/bin/phastaC.exe
  WORKING_DIRECTORY ${CDIR}
)
set(cmd 
  ${PHASTA_BINARY_DIR}/bin/checkphasta 
  ${CDIR}/2-procs_case-Posix/ 
  ${CDIR}/2-procs_case-Posix_ref/ 
  0 1e-6)
add_test(
  NAME compareCompressible-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME compressibleRestart-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${PHASTA_BINARY_DIR}/bin/phastaC.exe
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME compareCompressibleRestart-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
