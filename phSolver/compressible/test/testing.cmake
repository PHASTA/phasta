set(CDIR ${CASES}/compressible)
add_test(copyInpCfg 
  cp ${CMAKE_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
add_test(resetNumStart 
  cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
add_test(
  NAME compressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 2 ${CMAKE_BINARY_DIR}/bin/phastaC.exe
  WORKING_DIRECTORY ${CASES}/compressible
)
add_test(compareCompressible
  diff ${CDIR}/restart-dat.2.1.ref ${CDIR}/2-procs_case/restart-dat.2.1)
