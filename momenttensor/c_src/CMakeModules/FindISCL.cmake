# include(FindLibraryWithDebug)
#if (NOT ISCL_INCLUDE_DIR AND ISCL_LIBRARY)
  if (ISCL_INCLUDE_DIR AND ISCL_LIBRARY)
    set(ISCL_FIND_QUIETLY TRUE)
  endif (ISCL_INCLUDE_DIR AND ISCL_LIBRARY)
  find_path(ISCL_INCLUDE_DIR
    NAMES iscl_shared iscl_static
    HINTS /usr/include /usr/local/include $ENV{ISCL_DIR}/include
  )
  find_library(ISCL_LIBRARY
    NAMES iscl
    HINTS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 $ENV{ISCL_DIR}/lib
  )
  if (ISCL_INCLUDE_DIR)
    message("Found ISCL include")
  endif (ISCL_INCLUDE_DIR)
  if (ISCL_LIBRARY)
    message("Found ISCL library")
  endif (ISCL_LIBRARY)
  include(FindPackageHandleStandardArgs)
#endif()
find_package_handle_standard_args(ISCL DEFAULT_MSG
                                  ISCL_INCLUDE_DIR ISCL_LIBRARY)

mark_as_advanced(ISCL_INCLUDE_DIR ISCL_LIBRARY)
