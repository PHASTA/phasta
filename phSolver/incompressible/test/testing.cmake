set(CDIR ${CASES}/incompressible)
add_test(copyInpCfg
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
add_test(linkProcsDir-sync
  ln -snf ${CDIR}/4-procs_case-SyncIO-2 ${CDIR}/4-procs_case)
add_test(incompressibleResetNumStart-sync
  cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
add_test(
  NAME incompressible-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CASES}/incompressible
)
add_test(compareIncompressible-sync
  diff ${CDIR}/restart-dat.4.1.syncio.ref ${CDIR}/4-procs_case/restart-dat.4.1)
add_test(
  NAME incompressibleRestart-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CASES}/incompressible
)
add_test(compareIncompressibleRestart-sync
  diff ${CDIR}/restart-dat.8.1.syncio.ref ${CDIR}/4-procs_case/restart-dat.8.1)

add_test(linkProcsDir-posix
  ln -snf ${CDIR}/4-procs_case-Posix ${CDIR}/4-procs_case)
add_test(incompressibleResetNumStart-posix
  cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
add_test(
  NAME incompressible-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CASES}/incompressible
)
add_test(
  NAME incompressibleRestart-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CASES}/incompressible
)
