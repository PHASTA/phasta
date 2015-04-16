find_path(MESHSIM_INCLUDE_DIRS MeshSim.h HINTS "${MESHSIM_LIBRARY_DIR}/../../include" /usr/local/meshSim/latest /net/common/meshSim/latest)
find_path(MESHSIM_LIBRARY_DIR libSimMeshing.a HINTS /usr/local/meshSim/latest /net/common/meshSim/latest)
find_library(MESHSIM_LIBRARY libSimMeshing.a HINTS ${MESHSIM_LIBRARY_DIR} /usr/local/meshSim/latest /net/common/meshSim/latest)
find_library(PARASOLID_LIBRARY libpskernel.a HINTS ${MESHSIM_LIBRARY_DIR}/psKrnl /usr/local/meshSim/latest /net/common/meshSim/latest)

mark_as_advanced(PARASOLID_LIBRARY)
mark_as_advanced(MESHSIM_LIBRARY)

set(MESHSIM_LIBRARIES ${MESHSIM_LIBRARY}) #TODO: get the others too
