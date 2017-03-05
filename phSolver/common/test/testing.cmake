macro(common_parallel_test name procs dir exe)
  set(tname common_${name})
  add_test(
    NAME ${tname}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${procs} ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_common")
endmacro(common_parallel_test)

common_parallel_test(readHeader 4 ${CASES}/incompressible
  ${PHASTA_BINARY_DIR}/bin/phIOreadheader 2)
if(PHASTA_CHEF_ENABLED)
  common_parallel_test(writeReadZeroSz 2 ${CASES}/incompressible
    ${PHASTA_BINARY_DIR}/bin/phIOwriteReadZeroSz  1)
  common_parallel_test(writeFields 1 ${CASES}/incompressible
    ${PHASTA_BINARY_DIR}/bin/phIOwriteFields  1)
endif()
common_parallel_test(readIlwork 4
  ${CASES}/crossflow/4-1chef/4-procs_case
  ${PHASTA_BINARY_DIR}/bin/phIOreadIlwork
  . 0 0 foo)
common_parallel_test(readHeaderMultiTopo 4
  ${CASES}/crossflow/4-1chef/4-procs_case
  ${PHASTA_BINARY_DIR}/bin/phIOposixMultiTopo)
common_parallel_test(readDatablock 4 ${CASES}/incompressible
  ${PHASTA_BINARY_DIR}/bin/phIOreaddatablock 2)
common_parallel_test(write 4 ${CASES}
  ${PHASTA_BINARY_DIR}/bin/phIOwrite 2)
common_parallel_test(readFtn 4 ${CASES}/incompressible/
  ${PHASTA_BINARY_DIR}/bin/phIOreadFtn)
if(HAS_VALGRIND)
  common_parallel_test(readFtnVG 4 ${CASES}/incompressible/
    valgrind --log-file=vg.%p --leak-check=yes
    ${PHASTA_BINARY_DIR}/bin/phIOreadFtn)
endif(HAS_VALGRIND)
common_parallel_test(writeFtn 4 ${CASES}
  ${PHASTA_BINARY_DIR}/bin/phIOwriteFtn)
if(HAS_VALGRIND)
  common_parallel_test(
    writeFtnVG 4 ${CASES}/incompressible
    valgrind --log-file=vg.%p --leak-check=yes
    ${PHASTA_BINARY_DIR}/bin/phIOwriteFtn)
endif(HAS_VALGRIND)
