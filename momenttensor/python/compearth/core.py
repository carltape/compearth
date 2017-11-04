#!/usr/bin/env python

import numpy as np
import obspy

from obspy.core.event.base import _event_type_class_factory
from obspy.core.event.header import ATTRIBUTE_HAS_ERRORS
from obspy.core.event.source import Tensor as TensorBase

from .change_basis import change_basis
from .convert import cmt2tt, cmt2tt15




class Tensor(TensorBase):
    """ 
    Adds experimental methods to obspy MomentTensor
    """

    def change_basis(self, basis_type):
        # choose from five different basis conventions
        return change_basis(_vec(self.tensor), 1, basis_type)


    def convert(self, parameter_type):
        # plugin system for user-supplied parameterizations
        try: 
            plugins = __import__('compearth.plugins.convert')
            parameter_map = getattr(plugins, 'parameter_type')
        except:
            raise ParameterError

        try:
            return parameter_map(self) 
        except:
            raise Exception

    def convert_2012(self):
        gamma, delta, M0, kappa, theta, sigma = cmt2tt(_vec(self))
        return Tensor2012(
            gamma=gamma,
            delta=delta,
            M0=M0,
            kappa=kappa,
            theta=theta,
            sigma=sigma)

    def convert_2015(self):
        rho, v, w, kappa, sigma, h = cmt2tt15(_vec(self))
        return Tensor2015(
            rho=rho,
            v=v,
            w=w,
            kappa=kappa,
            sigma=sigma,
            h=h)


Tensor2012 = _event_type_class_factory(
    "Tensor2012",
    class_attributes=[("gamma", float, ATTRIBUTE_HAS_ERRORS),
                      ("delta", float, ATTRIBUTE_HAS_ERRORS),
                      ("M0", float, ATTRIBUTE_HAS_ERRORS),
                      ("kappa", float, ATTRIBUTE_HAS_ERRORS),
                      ("theta", float, ATTRIBUTE_HAS_ERRORS),
                      ("sigma", float, ATTRIBUTE_HAS_ERRORS)])


Tensor2015 = _event_type_class_factory(
    "Tensor2015",
    class_attributes=[("rho", float, ATTRIBUTE_HAS_ERRORS),
                      ("v", float, ATTRIBUTE_HAS_ERRORS),
                      ("w", float, ATTRIBUTE_HAS_ERRORS),
                      ("kappa", float, ATTRIBUTE_HAS_ERRORS),
                      ("sigma", float, ATTRIBUTE_HAS_ERRORS),
                      ("h", float, ATTRIBUTE_HAS_ERRORS)])



def _vec(M):
    return np.array([
        M.m_rr,
        M.m_tt,
        M.m_pp,
        M.m_rt,
        M.m_rp,
        M.m_tp])


def _mat(M):
    m = _vec(M)
    return np.array(([[m[0], m[3], m[4]],
                      [m[3], m[1], m[5]],
                      [m[4], m[5], m[2]]]))

