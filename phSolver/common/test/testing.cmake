add_test(
  NAME readHeader
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phIOreadheader 2
  WORKING_DIRECTORY ${CASES}/incompressible
)
add_test(
  NAME readDatablock
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phIOreaddatablock 2
  WORKING_DIRECTORY ${CASES}/incompressible
)
add_test(
  NAME write
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phIOwrite 2
  WORKING_DIRECTORY ${CASES}
)
add_test(
  NAME writeFtn
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 valgrind --log-file=vg.%p ${PHASTA_BINARY_DIR}/bin/phIOwriteFtn
  WORKING_DIRECTORY ${CASES}/incompressible/4-procs_case-SyncIO-2
)