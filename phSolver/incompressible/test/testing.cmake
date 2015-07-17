set(CDIR ${CASES}/incompressible)
add_test(copyInpCfg
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
add_test(linkProcsDir
  ln -f -s ${CDIR}/4-procs_case-SyncIO-2 ${CDIR}/4-procs_case)
add_test(resetNumStart
  cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
add_test(
  NAME incompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CASES}/incompressible
)
add_test(compareIncompressible
  diff ${CDIR}/restart-dat.4.1.syncio.ref ${CDIR}/4-procs_case/restart-dat.4.1)
add_test(
  NAME incompressibleRestart
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CASES}/incompressible
)
add_test(compareIncompressibleRestart
  diff ${CDIR}/restart-dat.8.1.syncio.ref ${CDIR}/4-procs_case/restart-dat.8.1)
