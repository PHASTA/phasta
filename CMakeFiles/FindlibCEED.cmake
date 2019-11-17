find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
	pkg_check_modules(CEEDPKG "ceed")
else(PKG_CONFIG_FOUND)
	message(FATAL_ERROR "Using libCEED without pkg-config not yet implemented")
endif(PKG_CONFIG_FOUND)

if(CEEDPKG_STATIC_LIBRARY_DIRS)
	#Prefer static if we can
	#This isn't really supported and we'll probably get
	#shared anyway but it's worth a try
	set(CEED_LIBRARY_DIRS ${CEEDPKG_STATIC_LIBRARY_DIRS})
	set(CEED_INCLUDE_DIRS ${CEEDPKG_STATIC_INCLUDE_DIRS})
	set(CEED_LIBRARIES ${CEEDPKG_STATIC_LIBRARIES})
else()
	set(CEED_LIBRARY_DIRS ${CEEDPKG_LIBRARY_DIRS})
	set(CEED_INCLUDE_DIRS ${CEEDPKG_INCLUDE_DIRS})
	set(CEED_LIBRARIES ${CEEDPKG_LINK_LIBRARIES})
endif()

if(CEED_LIBRARIES)
	set(CEED_FOUND TRUE)
endif()
