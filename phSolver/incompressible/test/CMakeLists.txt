macro(ic_parallel_test name procs dir exe)
  set(tname incompressible_${solver}_${name})
  add_test(
    NAME ${tname}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${procs} ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_incompressible")
endmacro(ic_parallel_test)

macro(ic_serial_test name exe)
  set(tname incompressible_${solver}_${name})
  add_test( NAME ${tname} COMMAND ${exe} ${ARGN} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_incompressible")
endmacro(ic_serial_test)

set(CDIR ${CASES}/incompressible)

ic_serial_test(copyInpCfg
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

if(PHASTA_USE_SVLS)
  set(solver svls)
  ic_serial_test(solverInp ln -snf ${CDIR}/solver.inp.svls ${CDIR}/solver.inp)
  include(testing.cmake)
endif()

set(solver leslib)
ic_serial_test(solverInp ln -snf ${CDIR}/solver.inp.libles ${CDIR}/solver.inp)
include(testing.cmake)
