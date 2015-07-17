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
  NAME writeheader
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phIOwriteheader 2
  WORKING_DIRECTORY ${CASES}
)
