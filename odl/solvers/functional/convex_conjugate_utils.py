# Copyright 2014-2016 The ODL development group
#
# This file is part of ODL.
#
# ODL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ODL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ODL.  If not, see <http://www.gnu.org/licenses/>.

"""REWRITE!

Factory functions for creating proximal operators.

Functions with ``cconj`` mean the proximal of the convex conjugate and are
provided for convenience.

For more details see :ref:`proximal_operators` and references therein. For
more details on proximal operators including how to evaluate the proximal
operator of a variety of functions see [PB2014]_. """


# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

from odl.solvers import Functional
from odl import LinearSpaceVector

# __all__ = ('Functional',)


def convex_conjugate_translation(convex_conj_f, y):
    """Calculate the convex conjugate functional of the translated function
    F(x - y).
    
    This is calculated according to the rule

        (F( . - y))^* (x) = F^*(x) - <y, x>

    where ``y`` is the translation of the argument.

    Parameters
    ----------
    convex_conj_f : `Functional`
        Function corresponding to F^*.
    
    y : Element in domain of F^*.

    Returns
    -------
    ConvexConjugateTranslation : `Functional`
        Functional corresponding to (F( . - y))^*

    Notes
    -----
    For reference on the identity used, see [KP2015]_.
   
    """

    #TODO: Add the reference    
    
    if y is not None and not isinstance(y, LinearSpaceVector):
        raise TypeError('vector {!r} not None or a LinearSpaceVector instance.'
                        ''.format(y))
    
    if not y in convex_conj_f.domain:
        raise TypeError('vector {} not in the domain of the functional {}.'
                        ''.format(y, convex_conj_f.domain))
    
    class ConvexConjugateTranslation(Functional):
        def __init__(self, convex_conj_f, y):
            super().__init__(domain=convex_conj_f.domain,
                             linear=convex_conj_f.is_linear,
                             smooth=convex_conj_f.is_smooth,
                             concave=convex_conj_f.is_concave,
                             convex=convex_conj_f.is_convex)
            
            self.orig_convex_conj_f = convex_conj_f
            self.y = y
            
            # The Lipschitz constant for the gradient can be bounded, by using
            # triangle inequality. However: is it the tightest bound?

        def _call(self, x, out):
            out[:] = convex_conj_f(x) + x.inner(self.y)
            
            
        def gradient(self, x, out):
            out[:] = self.orig_convex_conj_f.gradient(x) + self.y

        #TODO: Add this when the proximal frame-work is added to the functional
#        def proximal(self, sigma=1.0):
#            """Return the proximal operator of the functional.
#    
#            Parameters
#            ----------
#            sigma : positive float, optional
#                Regularization parameter of the proximal operator  
#    
#            Returns
#            -------
#            out : Operator 
#                Domain and range equal to domain of functional 
#            """
#            raise NotImplementedError

        #TODO: Add this when convex conjugate of a linear perturbation has been
        # added
#        def conjugate_functional(self):
#            """Convex conjugate functional of the functional.
#    
#            Parameters
#            ----------
#            none
#            
#            Returns
#            -------
#            out : Functional 
#                Domain equal to domain of functional 
#            """
#            raise NotImplementedError


        def derivative(self, point):

            class DerivativeOperator(Functional):
                def __init__(self):
                    super().__init__(self.orig_convex_conj_f.domain,
                                     linear=True)
    
                self.point=point
    
                def _call(self, x):
                    return x.inner(self.y +
                                   self.orig_convex_conj_f.gradient(point)) 
            
            return DerivativeOperator()
            
    return ConvexConjugateTranslation(convex_conj_f, y)
    
    
    
    
    
    
    
    
    
    
def convex_conjugate_arg_scaling(convex_conj_f, scaling):
    """Calculate the convex conjugate of function F(x * scaling).

    -------------------- REWRITE BELOW -----------------------
    This is calculated according to the rule

        prox[s * F( . * scaling)](x) =
        1/scaling * prox[s * scaling^2 * F ](x * scaling)

    where ``scaling`` is the scaling parameter, and ``s`` is the step size.

    Parameters
    ----------
    prox_factory : `callable`
        A factory function that, when called with a step size, returns the
        proximal operator of ``F``
    scaling : `float`
        Scaling parameter

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    Notes
    -----
    For reference on the identity used, see [CP2011c]_.
    """

    scaling = float(scaling)
    if scaling == 0:
        return proximal_zero(prox_factory(1.0).domain)

    def arg_scaling_prox_factory(step_size):
        """Create proximal for the translation with a given step_size.

        Parameters
        ----------
        step_size : positive `float`
            Step size parameter

        Returns
        -------
        proximal : `Operator`
            The proximal operator of ``s * F( . * a)`` where ``s`` is the
            step size
        """
        prox = prox_factory(step_size * scaling ** 2)
        return (1 / scaling) * prox * scaling

    return arg_scaling_prox_factory    
    
