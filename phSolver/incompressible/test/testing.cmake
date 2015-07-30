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
  WORKING_DIRECTORY ${CDIR}
)
set(cmd 
  ${PHASTA_BINARY_DIR}/bin/checkphasta 
  ${CDIR}/4-procs_case-SyncIO-2/ 
  ${CDIR}/4-procs_case-SyncIO-2_ref/ 
  2 1e-6)
add_test(
  NAME compareIncompressible-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME incompressibleRestart-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME compareIncompressibleRestart-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)

add_test(linkProcsDir-posix
  ln -snf ${CDIR}/4-procs_case-Posix ${CDIR}/4-procs_case)
add_test(incompressibleResetNumStart-posix
  cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
add_test(
  NAME incompressible-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CDIR}
)
set(cmd 
  ${PHASTA_BINARY_DIR}/bin/checkphasta 
  ${CDIR}/4-procs_case-Posix/ 
  ${CDIR}/4-procs_case-Posix_ref/ 
  0 1e-6)
add_test(
  NAME compareIncompressible-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME incompressibleRestart-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME compareIncompressibleRestart-posix
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
