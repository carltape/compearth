#!/usr/bin/python3
"""
Python interface to the underlying C library.

Copyright: ISTI distributed under the MIT license.
"""
import os
from ctypes import cdll
from ctypes import c_int
from ctypes import c_double
from ctypes import byref
from ctypes import POINTER
from numpy import zeros
from numpy import array
from numpy import reshape
from numpy import float64
from numpy import ascontiguousarray 
from numpy import linspace
from numpy.random import rand

class compearth:
    def __init__(self,
                 compearth_path=os.environ['LD_LIBRARY_PATH'].split(os.pathsep),
                 compearth_library='libcompearth_shared.so'):
        lfound = False
        for path in compearth_path:
            celib = os.path.join(path, compearth_library)
            if (os.path.isfile(celib)):
                lfound = True
                break
        if (lfound):
            ce = cdll.LoadLibrary(celib) 
        else:
            print("Couldn't find libcompearth_shared")
        # Make the interfaces
        ce.compearth_h2theta.argtypes = (c_int,
                                         POINTER(c_double),
                                         POINTER(c_double))
        ce.compearth_theta2h.argtypes = (c_int,
                                         POINTER(c_double),
                                         POINTER(c_double))
        ce.compearth_gamma2v.argtypes = (c_int,
                                         POINTER(c_double),
                                         POINTER(c_double))
        ce.compearth_v2gamma.argtypes = (c_int,
                                         POINTER(c_double), 
                                         POINTER(c_double))
        self.ce = ce
        return 

    def __enter__(self):
        return self

    def __exit__(self):
        return

    def __arrayToFloat64Pointer__(self, array):
        array = ascontiguousarray(array, float64)
        arrayPtr = array.ctypes.data_as(POINTER(c_double))
        return arrayPtr 

    def __allocFloat64Pointer__(self, n):
        array = ascontiguousarray(zeros(n), float64)
        arrayPtr = array.ctypes.data_as(POINTER(c_double))
        return array, arrayPtr 
  
    def gamma2v(self, gamma):
        """
        Computes rectlinear coordinate v from the lune longitude with 
        Equation 24b of Tape and Tape, 2015.
 
        Input
        ----- 
        gamma : array_like
           Lune longitudes in radians such that gamma is in the range
           [-pi/6, pi/6].

        Output
        ------
        v : array_like
           Rectlinear coordinate v such that v is in the range [-1/3, 1/3].
        """
        n = len(gamma)
        gammaPtr = self.__arrayToFloat64Pointer__(gamma)
        v, vPtr = self.__allocFloat64Pointer__(n)
        self.ce.compearth_gamma2v(n, gammaPtr, vPtr)
        return v

    def v2gamma(self, v):
        """
        Computes the lune longitude from the point v.
 
        Input
        -----
        v : array_like
           v in rectilinear space s.t. v is in the range [-1/3, 1/3].

        Output
        ------
        gamma : array_like
           Lune longitude in radians s.t. gamma in the range [-pi/6, pi/6].
        """
        n = len(v)
        vPtr = self.__arrayToFloat64Pointer__(v)
        gamma, gammaPtr = self.__allocFloat64Pointer__(n)
        self.ce.compearth_v2gamma(n, vPtr, gammaPtr)
        return gamma 

if __name__ == "__main__":
   ce = compearth()
   v = (rand(55) - 0.5)*2/3.
   gamma = ce.v2gamma(v)
   vn = ce.gamma2v(gamma)
   print(v - vn)
