add_executable(phIOreadheader phIOreadheader.cc)
include_directories("${PHASTA_BINARY_DIR}")
target_link_libraries(phIOreadheader common phastaIO)
include(testing.cmake)
