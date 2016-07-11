macro(c_parallel_test name procs dir exe)
  set(tname ${case}_${solver}_${name})
  add_test(
    NAME ${tname}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${procs} ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_compressible")
endmacro(c_parallel_test)

macro(c_serial_test name exe)
  set(tname ${case}_${solver}_${name})
  add_test( NAME ${tname} COMMAND ${exe} ${ARGN} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_compressible")
endmacro(c_serial_test)

foreach(testName compressible compressibleShockTube)
  set(CDIR ${CASES}/${testName})
  set(case ${testName})
  
  c_serial_test(inpCfg cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
  
  if(PHASTA_USE_PETSC)
    set(solver petsc)
    c_serial_test(solverInp ln -snf ${CDIR}/solver.inp.petsc ${CDIR}/solver.inp)
    include(${testName}.cmake)
  endif()
  
  set(solver native)
  c_serial_test(solverInp ln -snf ${CDIR}/solver.inp.native ${CDIR}/solver.inp)
  include(${testName}.cmake)
endforeach()
