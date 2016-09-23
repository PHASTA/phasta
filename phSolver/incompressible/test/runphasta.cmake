macro(cmd dir exe)
  message("${exe} ${ARGN}")
  execute_process(
    COMMAND ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir} 
    OUTPUT_VARIABLE out
    ERROR_VARIABLE out
    RESULT_VARIABLE res
    )
  message("${out}")
  if(res)
    message(FATAL_ERROR "Error running ${exe}")
  else()
    message("Success")
  endif()
endmacro()

cmd(${WORKDIR} ln -snf ${CASEDIR} ${TGTCASEDIR})
if(DEFINED NUMSTART )
  cmd(${WORKDIR} cp ${NUMSTART} ${TGTCASEDIR}/numstart.dat)
endif()
cmd(${WORKDIR} ${MPIRUN} ${MPIRUN_PROCFLAG} ${NUMPROCS} ${EXE})
