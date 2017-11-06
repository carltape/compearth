#!/usr/bin/env python

import numpy as np

from compearth.convert import cmt2tt, cmt2tt15, tt2cmt, tt152cmt

EPSVAL = 1.e-6


# must be in up-south-east (GCMT) convention
M = np.array([
    1.006279239004, # m11
    0.737173428960, # m22
    0.558314768020, # m33
   -0.231591759935, # m12
   -0.111675288138, # m13
    0.004991096805, # m23
    ])


def test_2012():
    M1 = M

    # convert to 2012 parameters and back
    gamma, delta, M0, kappa, theta, sigma = cmt2tt(M1)
    M2 = tt2cmt(gamma, delta, M0, kappa, theta, sigma)

    e = np.linalg.norm(M1-M2)
    if e > EPSVAL:
        print '||M1 - M2|| = %e' % e
        raise Exception


def test_2015():
    M1 = M

    # convert to 2015 parameters and back
    rho, v, w, kappa, sigma, h = cmt2tt15(M1)
    M2 = tt152cmt(rho, v, w, kappa, sigma, h)

    e = np.linalg.norm(M1-M2)
    if e > EPSVAL:
        print '||M1 - M2|| = %e' % e
        raise Exception



if __name__ == '__main__':
    try:
        test_2012()
    except:
        print 'Test 1 of 2...FAILED'
    else:
        print 'Test 1 of 2...SUCCESS'

    try:
        test_2015()
    except:
        print 'Test 2 of 2...FAILED'
    else:
        print 'Test 2 of 2...SUCCESS'

