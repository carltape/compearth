# About

This is a C translation of some of Carl Tape's compearth Matlab scripts.

# Prerequisites

An effort was made to make this library largely free of many dependencies.  
However, 

* CMake
* A valid C11 compiler
* BLAS, CBLAS, LAPACK, and LAPACKE along with header files
* Alternatively, BLAS, CBLAS, LAPACK, LAPACKE, and the requisite header files are available in Intel's [MKL](https://software.intel.com/en-us/mkl).
  
The purpose of CMake is to generate valid Makefiles in a cross-platform way.  Unfortunately, I find it infuriating to use.  To mitigate some of the frustration I usually create a compilation script.  For example, to configure with MKL and the Intel compiler I would do something like:

    #!/bin/bash
    /usr/bin/cmake ./ -DCMAKE_BUILD_TYPE=DEBUG \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_C_FLAGS="-g3 -O2 -qopenmp -Wall -Wextra -Wcomment -Wcheck" \
    -DCOMPEARTH_USE_MKL=TRUE \
    -DMKL_LIBRARY="/opt/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64_lin/libmkl_sequential.so;/opt/intel/mkl/lib/intel64_lin/libmkl_core.so" \
    -DMKL_INCLUDE_DIR=/opt/intel/mkl/include



