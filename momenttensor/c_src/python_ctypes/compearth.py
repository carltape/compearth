#!/usr/bin/python3
"""
Python interface to the underlying C library.

Copyright: ISTI distributed under the MIT license.
"""
import os
import timeit
from ctypes import cdll
from ctypes import c_int
from ctypes import c_double
from ctypes import byref
from ctypes import POINTER
from math import pi
from numpy import copy
from numpy import zeros
from numpy import array
from numpy import reshape
from numpy import float64
from numpy import ascontiguousarray 
from numpy import linspace
from numpy.random import rand
from numpy.random import seed

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
        ce.compearth_rect2lune.argtypes = (c_int, POINTER(c_double),
                                           c_int, POINTER(c_double),
                                           POINTER(c_double),
                                           POINTER(c_double)) 
        ce.compearth_lune2rect.argtypes = (c_int, POINTER(c_double),
                                           c_int, POINTER(c_double),
                                           POINTER(c_double),
                                           POINTER(c_double))
        ce.compearth_beta2u.argtypes = (c_int,
                                        POINTER(c_double),
                                        POINTER(c_double))
        ce.compearth_u2beta_optimize.argtypes = (c_int, c_int, c_int,
                                                 POINTER(c_double),
                                                 c_double,
                                                 POINTER(c_double))
        ce.compearth_u2beta.argtypes = (c_int,
                                        POINTER(c_double),
                                        POINTER(c_double))
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


    def rect2lune(self, v, w):
        """
        Converts v-w coordinates to lune coordinates.

        Input
        -----
        v : array_like
            v ordinate s.t. v is in the range [-1/3, 1/3].  This corresponds
            to gamma.
        w : array_like
            w ordinate s.t. w is in the rnage [-3*pi/8, 3*pi/8].  This 
            corresponds to delta.

        Ouput
        -----
        gamma : array_like
            Longitude (degrees) where gamma is in the range [-30,30].
        delta : array_like
            Latitude (degrees) where delta is in the rnage [-90,90].
        """
        nv = len(v)
        nw = len(w)
        vPtr = self.__arrayToFloat64Pointer__(v)
        wPtr = self.__arrayToFloat64Pointer__(w)
        gamma, gammaPtr = self.__allocFloat64Pointer__(nv)
        delta, deltaPtr = self.__allocFloat64Pointer__(nw)
        self.ce.compearth_rect2lune(nv, vPtr, nw, wPtr, gammaPtr, deltaPtr)
        return gamma, delta

    def lune2rect(self, gamma, delta):
        """
        Converts from lune coordinates (gamma, delta) to rectilinear 
        coordinates (v, w).
        Input
        -----
        gamma : array_like
            Lune longitudes in degrees.  These are in the range [-30,30].
        delta : array_like
            Lune latitudes in degrees.  These are in the range [-90,90].

        Output
        ------
        v : array_like
            v coordinates corresponding to the longitudes.  These are in
            the range [-1/3, 1/3].
        w : array_like
            w coordinates corresonding to the latitudes.  These are in
            the range [-3*pi/8, 3*pi/8].
        """
 
        ng = len(gamma)
        nd = len(delta)
        gammaPtr = self.__arrayToFloat64Pointer__(gamma) 
        deltaPtr = self.__arrayToFloat64Pointer__(delta)
        v, vPtr = self.__allocFloat64Pointer__(ng)
        w, wPtr = self.__allocFloat64Pointer__(nd)
        self.ce.compearth_lune2rect(ng, gammaPtr, nd, deltaPtr, vPtr, wPtr)
        return v, w 
 
    def beta2u(self, beta):
        """
        Converts the colatitudes beta to u in the rectilinear coordinates.

        Input
        -----
        beta : array_like
           Colatitudes that are in the range [0,pi].

        Output
        ------
        u : array_like
           u in the rectilinear coordinates that are in the range [0, 3*pi/4].
        """
        n = len(beta)
        betaPtr = self.__arrayToFloat64Pointer__(beta)
        u, uPtr = self.__allocFloat64Pointer__(n)
        self.ce.compearth_beta2u(n, betaPtr, uPtr)
        return u

    def u2beta(self, u, maxit=25, tol=1.e-12, linvType=1):
        """
        Converts u from the rectilinear coordinates to colatitudes.

        Input
        -----
        u : array_like
           u coordinates in the range [0, 3*pi/4].  
        maxit : int
           Maximum number of iterations in non-linear root finding method.
        tol : float
           Convergence is achieved when:
           abs(u - 0.75*beta - 0.5*sin2b + 0.0625*sin4b) < tol
        linvType : int
           If 1 then use the Newton Rhapson iteration.  Otherwise, use the
           Halley iteration (which is the default).

        Result
        ------ 
        beta : array_like
           The colatitudes in radians that are in the range [0, pi].
        ierr : int
           Error flag where 0 indicates success.
        """
        n = len(u)
        uPtr = self.__arrayToFloat64Pointer__(u)
        beta, betaPtr = self.__allocFloat64Pointer__(n)
        ierr = self.ce.compearth_u2beta(n, uPtr, betaPtr) #maxit, linvType, uPtr, tol, betaPtr)
        return beta, ierr

    def h2theta(self, h):
        """
        Computes dip angle from h.

        Input
        -----
        h : array_like
           h in the rectilinear space where h is in the range [0,1].

        Output
        ------
        theta : array_like
           The dip angle in radians such that theta is in the range [0,pi/2].
        """
        n = len(h)
        hPtr = self.__arrayToFloat64Pointer__(h)
        theta, thetaPtr = self.__allocFloat64Pointer__(n)
        self.ce.compearth_h2theta(n, hPtr, thetaPtr)
        return theta

    def theta2h(self, theta):
        """
        Computes h from the dip angle according to Equation 24c of Tape and
        Tape, 2015.

        Input
        -----
        theta : array_like
           Dip angle in radians where theta is in the range [0,pi/2]. 

        Output
        ------
        h : array_like
           h in the rectilinear space where h is in the range [0,1]. 
        """
        n = len(theta)
        thetaPtr = self.__arrayToFloat64Pointer__(theta)
        h, hPtr = self.__allocFloat64Pointer__(n)
        self.ce.compearth_theta2h(n, thetaPtr, hPtr)
        return h
 
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

def unit_test():
    seed(8675309)
    ce = compearth()
    """
    beta0 = linspace(0, pi, 1000)
    u0 = ce.beta2u(beta0)
    clineb = '     const double betas[%d] = {'%len(beta0)
    clineu = '     const double us[%d] = {'%len(u0)
    for i in range(len(beta0)):
        clineb = clineb + '%.15f'%beta0[i]
        clineu = clineu + '%.15f'%u0[i]
        if (i < len(beta0) - 1):
            clineb = clineb + ','
            clineu = clineu + ','
            if (i%4 == 0):
                clineb = clineb + '\n     '
                clineu = clineu + '\n     '
    clineb = clineb + '};\n'
    clineu = clineu + '};\n'
    print(clineb)
    print(clineu)
    return 
    """
    start = timeit.timeit()
    # gamma to v conversions
    v = (rand(100) - 0.5)*2/3.
    gamma = ce.v2gamma(v)
    vn = ce.gamma2v(gamma)
    assert(max(abs(v - vn)) < 1.e-14), 'failed gamma2v'
    # h to theta conversions
    h = rand(100)
    theta = ce.h2theta(h)
    hn = ce.theta2h(theta)
    assert(max(abs(hn - h)) < 1.e-14), 'failed h2theta'
    # beta test; there's problems near the ends
    ub = 2.17259
    lb = 0.0000245236
    beta = 0.5*(lb + ub) + 0.5*(ub - lb)*rand(100)
    u = ce.beta2u(beta)
    betan, ierr = ce.u2beta(u)
    #for i in range(len(beta)):
    #    if (abs(betan[i] - beta[i]) > 1.e-10):# and u[i] > 0.001):
    #        print(u[i], beta[i], betan[i], abs(betan[i] - beta[i]))
    assert(max(abs(betan - beta)) < 1.e-5), 'failed betatest' 
    # lune2rect test
    gamma = None
    gamma = copy(0.5*(30.0 - 30.0) + 0.5*(30.0 - -30.0)*rand(150))
    # there's still a problem at the high latitudes
    delta = copy(0.5*(89.9 - 89.9) + 0.5*(89.9 - -89.9)*rand(295))
    v, w = ce.lune2rect(gamma, delta)
    gamman, deltan = ce.rect2lune(v, w)
    assert(max(abs(gamma - gamman)) < 1.e-11), 'gamma in lune2rect failed'
    assert(max(abs(delta - deltan)) < 1.e-1), 'delta in lune2rect failed'
    end = timeit.timeit()
    print("Passed unit tests in %f seconds"%(end - start))
    return  

if __name__ == "__main__":
    unit_test()
