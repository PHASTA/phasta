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
if(HAS_VALGRIND)
  set(vgcmd
    valgrind
    --log-file=vg.%p
    --leak-check=yes
    ${PHASTA_BINARY_DIR}/bin/phIOreadFtn
  )
  add_test(
    NAME readFtnVG
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${vgcmd}
    WORKING_DIRECTORY ${CASES}/incompressible/
  )
endif(HAS_VALGRIND)
add_test(
  NAME writeFtn
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${PHASTA_BINARY_DIR}/bin/phIOwriteFtn
  WORKING_DIRECTORY ${CASES}/
)
if(HAS_VALGRIND)
  set(vgcmd
    valgrind
    --log-file=vg.%p
    --leak-check=yes
    ${PHASTA_BINARY_DIR}/bin/phIOwriteFtn
  )
  add_test(
    NAME writeFtnVG
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${vgcmd}
    WORKING_DIRECTORY ${CASES}/incompressible/
  )
endif(HAS_VALGRIND)
