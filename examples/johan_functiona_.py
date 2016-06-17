# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:36:12 2016

@author: johan79
"""

from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

import numpy as np
import odl

# Discretization parameters
n = 4

# Discretized spaces
space = odl.uniform_discr([0, 0, 0], [1, 1, 1], [n, n, n])


print(space)

class L1Norm(odl.solvers.functional.Functional):
    def __init__(self, domain):
        super().__init__(domain=domain, linear=False, convex=True, concave=False, smooth=False, grad_lipschitz=np.inf)

    def _call(self, x):
        return np.abs(x).inner(self.domain.one())

    @property
    def gradient(x):
        raise NotImplementedError

    def proximal(self, sigma=1.0):
        functional = self

        class L1Proximal(odl.operator.Operator):
            def __init__(self):
                super().__init__(functional.domain, functional.domain,
                                 linear=False)
                self.sigma = sigma

            #TODO: Check that this works for complex x
            def _call(self, x):
                return np.maximum(np.abs(x)-sigma,0)*np.sign(x)

        return L1Proximal()
    @property
    def conjugate_functional(self):
        functional = self

        class L1Conjugate_functional(odl.solvers.functional.Functional):
            def __init__(self):
                super().__init__(functional.domain, linear=False)

            def _call(self, x):
                if np.max(np.abs(x)) > 1:
                    return np.inf
                else:
                    return 0

        return L1Conjugate_functional()

class L2Norm(odl.solvers.functional.Functional):
    def __init__(self, domain):
        super().__init__(domain=domain, linear=False, convex=True, concave=False, smooth=False, grad_lipschitz=np.inf)

    def _call(self, x):
        return np.sqrt(np.abs(x).inner(np.abs(x)))

    @property
    def gradient(self):
        functional = self

        class L2Gradient(odl.operator.Operator):
            def __init__(self):
                super().__init__(functional.domain, functional.domain,
                                 linear=False)

            def _call(self, x):
                return x/x.norm()

        return L2Gradient()

    def proximal(self, sigma=1.0):
        functional = self

        class L2Proximal(odl.operator.Operator):
            def __init__(self):
                super().__init__(functional.domain, functional.domain,
                                 linear=False)
                self.sigma = sigma

            #TODO: Check that this works for complex x
            def _call(self, x):
                return np.maximum(x.norm()-sigma,0)*(x/x.norm())

        return L2Proximal()

    @property
    def conjugate_functional(self):
        functional = self

        class L2Conjugate_functional(odl.solvers.functional.Functional):
            def __init__(self):
                super().__init__(functional.domain, linear=False)

            def _call(self, x):
                if x.norm() > 1:
                    return np.inf
                else:
                    return 0

        return L2Conjugate_functional()

class L2NormSquare(odl.solvers.functional.Functional):
    def __init__(self, domain):
        super().__init__(domain=domain, linear=False, convex=True, concave=False, smooth=True, grad_lipschitz=2)

    def _call(self, x):
        return np.abs(x).inner(np.abs(x))

    @property
    def gradient(self):
        functional = self

        class L2SquareGradient(odl.operator.Operator):
            def __init__(self):
                super().__init__(functional.domain, functional.domain,
                                 linear=False)

            def _call(self, x):
                return 2*x

        return L2SquareGradient()

    def proximal(self, sigma=1.0):
        functional = self

        class L2SquareProximal(odl.operator.Operator):
            def __init__(self):
                super().__init__(functional.domain, functional.domain,
                                 linear=False)
                self.sigma = sigma

            #TODO: Check that this works for complex x
            def _call(self, x):
                return x/3

        return L2SquareProximal()

    @property
    def conjugate_functional(self):
        functional = self

        class L2SquareConjugateFunctional(odl.solvers.functional.Functional):
            def __init__(self):
                super().__init__(functional.domain, linear=False)

            def _call(self, x):
                return x.norm()/4

        return L2SquareConjugateFunctional()





l1func = L1Norm(space)
l1prox = l1func.proximal(sigma=1.5)
l1conjFun = l1func.conjugate_functional


# Create phantom
phantom = odl.util.shepp_logan(space, modified=True)*5+1


onevector=space.one()*5

prox_phantom=l1prox(phantom)
l1conjFun_phantom = l1conjFun(phantom)

l2func=L2Norm(space)
l2prox = l2func.proximal(sigma=1.5)
l2conjFun = l2func.conjugate_functional
l2conjGrad = l2func.gradient

prox2_phantom=l2prox(phantom*10)
l2conjFun_phantom = l2conjFun(phantom/10)

l22=L2NormSquare(space)
prox22=l22.proximal(1)(phantom)

l22(phantom)
cf22=l22.conjugate_functional(phantom)

