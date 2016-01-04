
find_path(FMDB_INCLUDE_DIR FMDB.h HINTS "${FMDB_DIR}/include")
find_library(FMDB_LIBRARY libFMDB.a NAMES libFMDB.so HINTS "${FMDB_DIR}/lib")
find_path(FMDB_DIR include/FMDB.h NO_DEFAULT_PATH)
mark_as_advanced(FMDB_INCLUDE_DIR)
mark_as_advanced(FMDB_LIBRARY)

set(FMDB_LIBRARIES ${FMDB_LIBRARY}) #TODO: get the others too
