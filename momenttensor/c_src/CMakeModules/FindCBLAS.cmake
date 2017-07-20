# include(FindLibraryWithDebug)
if (CBLAS_INCLUDE_DIR AND CBLAS_LIBRARY AND BLAS_LIBRARY)
  set(CBLAS_FIND_QUIETLY TRUE)
endif (CBLAS_INCLUDE_DIR AND CBLAS_LIBRARY AND BLAS_LIBRARY)
find_path(CBLAS_INCLUDE_DIR
  NAMES cblas.h
  HINTS /usr/include $ENV{CBLASDIR}/include ${INCLUDE_INSTALL_DIR}
)
find_library(CBLAS_LIBRARY
  NAMES cblas
  HINTS /usr/lib /usr/lib64 $ENV{CBLASDIR}/lib
)
# Look for openBLAS first
find_library(BLAS_LIBRARY
  NAME openblas
  HINTS /usr/lib /usr/lib64 $ENV{BLASDIR}/lib
)

find_library(BLAS_LIBRARY
  NAMES blas
  HINTS /usr/lib /usr/lib64 $ENV{BLASDIR}/lib
)

#find_file(CBLAS_LIBRARY
#  libcblas.so.3
#  libcblas.a
#  libblas.so.3
#  libblas.so
#  libblas.a
#  PATHS
#  /usr/lib
#  $ENV{CBLASDIR}/lib
#  ${LIB_INSTALL_DIR}
#)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS DEFAULT_MSG
                                  CBLAS_INCLUDE_DIR CBLAS_LIBRARY BLAS_LIBRARY)
mark_as_advanced(CBLAS_INCLUDE_DIR CBLAS_LIBRARY BLAS_LIBRARY)
