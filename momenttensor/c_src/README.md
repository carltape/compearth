# About

This is a C translation of some of Carl Tape's compearth Matlab scripts.  These functions are particularly useful when discretizing uniform moment tensor spaces or performing moment tensor decompositions in a compiled language such as C or Fortran.

# Prerequisites

An effort was made to make this library largely free of many dependencies.  
However, 

* CMake
* A valid C11 compiler
* BLAS, CBLAS, LAPACK, and LAPACKE along with header files
* Alternatively, BLAS, CBLAS, LAPACK, LAPACKE, and the requisite header files are available at no cost in Intel's [MKL](https://software.intel.com/en-us/mkl).
* Doxygen.  This is only required if you wish to generate the API documentation.
  
The purpose of CMake is to generate valid Makefiles in a cross-platform way.  Unfortunately, I find it infuriating to use.  To mitigate some of the frustration I usually create a compilation script.  For example, to configure with MKL and the Intel compiler I would do something like:

    #!/bin/bash
    /usr/bin/cmake ./ -DCMAKE_BUILD_TYPE=DEBUG \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_C_FLAGS="-g3 -O2 -qopenmp -Wall -Wextra -Wcomment -Wcheck" \
    -DCOMPEARTH_USE_MKL=TRUE \
    -DMKL_LIBRARY="/opt/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64_lin/libmkl_sequential.so;/opt/intel/mkl/lib/intel64_lin/libmkl_core.so" \
    -DMKL_INCLUDE_DIR=/opt/intel/mkl/include

To generate the Doxygen documentation one can do the following

    doxygen Doxyfile

from the source root directory.

# TODO List

The library is still limited in some aspects.  The following  

* Generate/upload F90 interface modules
* Extend unit tests to chunks of matrices
* Compute intelligent chunk sizes

