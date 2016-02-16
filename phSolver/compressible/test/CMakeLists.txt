macro(c_parallel_test name procs dir exe)
  set(tname compressible_${solver}_${name})
  add_test(
    NAME ${tname}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${procs} ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_compressible")
endmacro(c_parallel_test)

macro(c_serial_test name exe)
  set(tname compressible_${solver}_${name})
  add_test( NAME ${tname} COMMAND ${exe} ${ARGN} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_compressible")
endmacro(c_serial_test)

set(CDIR ${CASES}/compressible)

c_serial_test(inpCfg cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

if(PHASTA_USE_PETSC)
  c_serial_test(solverInp ln -snf ${CDIR}/solver.inp.petsc ${CDIR}/solver.inp)
  set(solver petsc)
  include(testing.cmake)
endif()

c_serial_test(solverInp ln -snf ${CDIR}/solver.inp.native ${CDIR}/solver.inp)
set(solver native)
include(testing.cmake)
