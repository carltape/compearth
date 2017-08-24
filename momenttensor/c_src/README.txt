[About]

This is a C translation of some of Carl Tape's compearth Matlab scripts.

[Prerequisites]

(1) LAPACKE and CBLAS.  These are freely available in Intel's MKL.

The build relies on CMake.  I usually use a configuration script instead of
arguing with CMake.

For example, I'll build with MKL.

#!/bin/bash
/usr/bin/cmake ./ -DCMAKE_BUILD_TYPE=DEBUG \
-DCMAKE_INSTALL_PREFIX=./ \
-DCMAKE_C_COMPILER=icc \
-DCMAKE_C_FLAGS="-g3 -O2 -qopenmp -Wall -Wextra -Wcomment -Wcheck" \
-DCOMPEARTH_USE_MKL=TRUE \
-DMKL_LIBRARY="/opt/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64_lin/libmkl_sequential.so;/opt/intel/mkl/lib/intel64_lin/libmkl_core.so" \
-DMKL_INCLUDE_DIR=/opt/intel/mkl/include


