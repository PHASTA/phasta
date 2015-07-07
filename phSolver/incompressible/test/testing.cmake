set(CDIR ${CASES}/crossflow)
add_test(copyInpCfg 
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
add_test(resetNumStart 
  cp ${CDIR}/numstart.dat ${CDIR}/4-procs_case/numstart.dat)
add_test(
  NAME crossflow
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  WORKING_DIRECTORY ${CASES}/crossflow
)
add_test(compareCrossflow 
  diff ${CDIR}/restart-dat.4.1.ref ${CDIR}/4-procs_case/restart-dat.4.1)
