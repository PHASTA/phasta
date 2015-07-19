macro(setup_exe exename srcname)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} common phastaIO)
endmacro(setup_exe)

include_directories("${PHASTA_BINARY_DIR}")

setup_exe(phIOreadheader phIOreadheader.cc)
setup_exe(phIOreaddatablock phIOreaddatablock.cc)
setup_exe(phIOwrite phIOwrite.cc)
setup_exe(phIOwriteFtn phIOwrite.f)
set_property(TARGET phIOwriteFtn PROPERTY LINKER_LANGUAGE Fortran)
include(testing.cmake)
