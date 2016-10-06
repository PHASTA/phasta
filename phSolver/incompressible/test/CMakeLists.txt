set(CDIR ${CASES}/incompressible)

add_test(NAME incompressible_copyInpCfg
  COMMAND cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

if(PHASTA_USE_SVLS)
  set(solver svls)
  add_test(NAME incompressible_${solver}_solverInp
    COMMAND ln -snf ${CDIR}/solver.inp.svls ${CDIR}/solver.inp)
  include(testing.cmake)
endif()

set(solver leslib)
add_test(NAME incompressible_${solver}_solverInp
  COMMAND ln -snf ${CDIR}/solver.inp.libles ${CDIR}/solver.inp)
include(testing.cmake)
