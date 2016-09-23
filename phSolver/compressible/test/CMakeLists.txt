foreach(testName compressible compressibleShockTube)
  set(CDIR ${CASES}/${testName})
  set(case ${testName})

  add_test(NAME ${case}_inpCfg
    COMMAND cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

  if(PHASTA_USE_PETSC)
    set(solver petsc)
    add_test(NAME ${case}_${solver}_solverInp
      COMMAND ln -snf ${CDIR}/solver.inp.petsc ${CDIR}/solver.inp)
    include(${testName}.cmake)
  endif()

  set(solver native)
  add_test(NAME ${case}_${solver}_solverInp
    COMMAND ln -snf ${CDIR}/solver.inp.native ${CDIR}/solver.inp)
  include(${testName}.cmake)
endforeach()
