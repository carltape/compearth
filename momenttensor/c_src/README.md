# About

This is a C translation of some of Carl Tape's compearth Matlab scripts.  These functions are particularly useful when discretizing uniform moment tensor spaces or performing moment tensor decompositions in a compiled language such as C or Fortran.

# Prerequisites

An effort was made to make this library largely free of many dependencies.  However, there still remain some required items which are detailed in the following list

* [CMake](https://cmake.org/) v2.6 or greater.
* A valid C11 compiler
* [LAPACK(E)](http://www.netlib.org/lapack/) and [(C)BLAS](http://www.netlib.org/blas/) along with header files.  Vendor implementations are recommended when available.
* When using Intel processors BLAS, CBLAS, LAPACK, LAPACKE, and the requisite header files are available at Intel's [MKL](https://software.intel.com/en-us/mkl).  This library is available at no cost.
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/).  This is only required if you wish to generate the API documentation.

# Build Instructions
 
The purpose of CMake is to serve as a cross-platform Makefile generator.  CMake, however, can be difficult to use.  I usually find the most expedient strategy is to create configuration scripts.  For eample, to configure with clang one could run the following script in the root source directory 

## Configuring With System LAPACK(e) and (C)BLAS

    #!/bin/bash
    /usr/bin/cmake ./ -DCMAKE_BUILD_TYPE=DEBUG \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_C_FLAGS="-g3 -O2 -fopenmp -Weverything" \
    -DLAPACKE_INCLUDE_DIR=/usr/include \
    -DLAPACKE_LIBRARY=/usr/lib/liblapacke.so \
    -DLAPACK_LIBRARY=/usr/lib/liblapack.so \
    -DCBLAS_INCLUDE_DIR=/usr/include \
    -DCBLAS_LIBRARY=/usr/lib/libblas.so \
    -DBLAS_LIBRARY=/usr/lib/libblas.so

## Configuring With MKL

Alternatively, to configure with MKL and the Intel C and Fortran compilers (which are not required by MKL) I would do something like:

    #!/bin/bash
    /usr/bin/cmake ./ -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DCOMPEARTH_MAKE_FINTER=TRUE \
    -DCMAKE_C_FLAGS="-g3 -xHOST -O3 -ipo -qopenmp -Wall -Wextra -Wcomment -Wcheck" \
    -DCMAKE_Fortran_FLAGS="-g2 -O2" \
    -DCOMPEARTH_USE_MKL=TRUE \
    -DMKL_LIBRARY="/opt/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64_lin/libmkl_sequential.so;/opt/intel/mkl/lib/intel64_lin/libmkl_core.so" \
    -DMKL_INCLUDE_DIR=/opt/intel/mkl/include

## Python Interfacing With CTypes

With [ctypes](https://docs.python.org/3/library/ctypes.html) and [NumPy](http://www.numpy.org/) it is possible to wrap some components of the library for use directly from Python.  While this won't yield the fastest interfaces it will be able to directly call the underlying C library and produce a portable solution.  For this to work one must be sure to only link to shared libraries.  The magic happens for Linux users by scraping your LD_LIBRARY_PATH; thus for this to work it is important that libcompearth_shared.so be one of your LD_LIBRARY_PATH's.  This environment variable is configurable in your .bashrc or .cshrc file.

Additionally, if using MKL, it may be necessary to link to libmkl_avx2, libmkl_mc3, and libmkl_def.  Unfortnately, this will only be obvious should you see something like: Intel MKL FATAL ERROR: Cannot load libmkl_avx2.so or libmkl_def.so.
 
## Experimental Python Interfacing - this is going to be eliminated

With [SWIG](http://www.swig.org/) and [NumPy](http://www.numpy.org/) it may be possible to wrap some components of the library for use from Python.  In this instance a configuration script may look like

    #!/bin/bash
    /usr/bin/cmake ./ -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=icc \
    -DCOMPEARTH_BUILD_PYTHON_MODULE=TRUE \
    -DSWIG_EXECUTABLE=/usr/bin/swig \
    -DCMAKE_C_FLAGS="-g3 -xHOST -O3 -ipo -qopenmp -Wall -Wextra -Wcomment -Wcheck" \
    -DCOMPEARTH_USE_MKL=TRUE \
    -DMKL_LIBRARY="/opt/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64_lin/libmkl_sequential.so;/opt/intel/mkl/lib/intel64_lin/libmkl_core.so;/opt/intel/mkl/lib/intel64_lin/libmkl_mc3.so;/opt/intel/lib/intel64_lin/libmkl_def.so" \
    -DMKL_INCLUDE_DIR=/opt/intel/mkl/include

Note, it may be necessary to include libmkl_avx2, libmkl_mc3, and libmkl_def as needed to the MKL_LIBRARY field.  In the instance of a failure one may see a terse message from Python similar to:

Intel MKL FATAL ERROR: Cannot load libmkl_avx2.so or libmkl_def.so.

Unfortunately, this will make building the Python module and testing the Python module an iterative process.

# Building

After configuration the application can be built in the source root directory by typing

    make

and, assuming root privelages are not required, installed by typing 

    make install

The library can be tested by typing

    make test

The Doxygen documentation can be generated by typing

    doxygen Doxyfile

from the source root directory.

# TODO List

The library is still limited in some aspects.  The following  

* Add more interface modules.  Deal with NULLs.  Might need Fortran pointers.  Yick.
* Extend unit tests to chunks of matrices
* Compute intelligent chunk sizes
* Wrap more functions in ctypes
* Make a setup.py for the ctypes wrappers

# Credits

This was written by Ben Baker of [ISTI](http://www.isti.com/) and is a derivative of Carl Tape's compearth matlab library.

