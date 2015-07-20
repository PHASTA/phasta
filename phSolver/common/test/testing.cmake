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
  NAME readFtn
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phIOreadFtn
  WORKING_DIRECTORY ${CASES}/incompressible/
)
add_test(
  NAME readFtnVG
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 valgrind --log-file=vg.%p --leak-check=yes ${PHASTA_BINARY_DIR}/bin/phIOreadFtn
  WORKING_DIRECTORY ${CASES}/incompressible/
)
