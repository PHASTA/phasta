set(testLabel "phsolver_incompressible_${solver}")
add_test(NAME ${testLabel}_sync
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-SyncIO-2
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DNUMSTART=${CDIR}/numstart.dat
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )

add_test(NAME ${testLabel}_restart-sync
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-SyncIO-2
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )

add_test(NAME ${testLabel}_posix
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-Posix
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DNUMSTART=${CDIR}/numstart.dat
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )

add_test(NAME ${testLabel}_restart-posix
  COMMAND ${CMAKE_COMMAND}
  -DNAME=${casename}
  -DWORKDIR=${CDIR}
  -DCASEDIR=${CDIR}/4-procs_case-Posix
  -DTGTCASEDIR=${CDIR}/4-procs_case
  -DMPIRUN=${MPIRUN}
  -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
  -DEXE=${PHASTA_BINARY_DIR}/bin/phastaIC.exe
  -DNUMPROCS=4
  -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
  )
