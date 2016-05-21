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
# along with ODL. If not, see <http://www.gnu.org/licenses/>.

"""Factory functions for creating proximal operators.

For more details see :ref:`proximal_operators` and references therein. For
more details on proximal operators including how to evaluate the proximal
operator of a variety of functions see [PB2014]_. """

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

import numpy as np
import scipy as sp

from odl.discr.tensor_ops import PointwiseNorm
from odl.operator.operator import Operator
from odl.operator.default_ops import IdentityOperator
from odl.operator.pspace_ops import ProductSpaceOperator
from odl.space.pspace import ProductSpace


__all__ = ('combine_proximals', 'proximal_convexconjugate',
           'proximal_composition', 'proximal_zero',
           'proximal_clamp', 'proximal_nonnegativity',
           'proximal_l1', 'proximal_convexconjugate_l1',
           'proximal_l2', 'proximal_convexconjugate_l2',
           'proximal_l2_squared', 'proximal_convexconjugate_l2_squared',
           'proximal_convexconjugate_kl', 'proximal_variable_lp')


# TODO: remove diagonal op once available on master
def combine_proximals(factory_list):
    """Combine proximal operators into a diagonal product space operator.

    This assumes the functional to be separable across variables in order to
    make use of the separable sum property of proximal operators.

        prox[tau * (f(x) + g(y))](x, y) =
            (prox[tau * f](x), prox[tau * g](y))

    Parameters
    ----------
    factory_list : list of `callable`
        A list containing proximal operator factories

    Returns
    -------
    diag_op : `callable`
        Returns a diagonal product space operator factory to be initialized
        with the same step size parameter
    """

    def diagonal_operator(operators, dom=None, ran=None):
        """Broadcast argument to set of operators.

        Parameters
        ----------
        operators : array-like
            An array of `Operator`s
        dom : `ProductSpace`, optional
            Domain of the operator. If not provided, it is tried to be
            inferred from the operators. This requires each **column**
            to contain at least one operator.
        ran : `ProductSpace`, optional
            Range of the operator. If not provided, it is tried to be
            inferred from the operators. This requires each **row**
            to contain at least one operator.
        """

        indices = [range(len(operators)), range(len(operators))]
        shape = (len(operators), len(operators))
        op_matrix = sp.sparse.coo_matrix((operators, indices), shape)

        return ProductSpaceOperator(op_matrix, dom=dom, ran=ran)

    def make_diag(step_size):
        """Diagonal matrix of operators

        Parameters
        ----------
        step_size : positive `float`
            Step size parameter

        Returns
        -------
        diag_op : `Operator`
        """
        return diagonal_operator(
            [factory(step_size) for factory in factory_list])

    return make_diag


def proximal_convexconjugate(proximal):
    """ Calculate the proximal of the dual using Moreau decomposition

    The Moreau identity states that for any convex function ``F`` with
    convex conjugate ``F^*``, the proximals satisfy

        prox[a * F^*](x) + a * prox[F / a](x / a) = x

    where ``a`` is a scalar step size. Using this, we find the proximal of the
    convex conjugate

        prox[a * F^*](x) = x - a * prox[F / a](x / a)

    Note that since ``(F^*)^* = F``, this can be used to get the primal from
    the dual.

    Parameters
    ----------
    proximal : `callable`
        A factory function that when called with a stepsize returns the
        proximal operator of ``F``.

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized
    """

    def make_convexconjugate_proximal(step_size):
        """Create proximal for the dual with a given step_size

        Parameters
        ----------
        step_size : positive `float`
            Step size parameter

        Returns
        -------
        proximal : `Operator`
            The proximal operator of ``prox[a * f^*](x)``
        """
        prox_other = step_size * proximal(1.0 / step_size) * (1.0 / step_size)
        return IdentityOperator(prox_other.domain) - prox_other

    return make_convexconjugate_proximal


def proximal_composition(proximal, operator, mu):
    """ Calculate the proximal of functional a composed with unitary operator

    Given an linear `Operator` ``L`` with the property that for a scalar ``mu``

        L^*(L(x)) = mu * x

    And a convex function ``f``, The following identity holds

        prox[f * L](x) = x + 1/mu L^*(prox[mu * f](Lx) - Lx)

    This factory function implements this functionality.

    There is no simple formula for more general operators.

    Parameters
    ----------
    proximal : `callable`
        A factory function that when called with a stepsize returns the
        proximal operator of ``F``.
    operator : `Operator`
        The operator to compose the functional with
    mu : `Operator.field` element
        Scalar such that ``operator.adjoint * operator = mu``

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    Notes
    -----
    The function cannot verify that the identity holds, the user needs to
    verify this.
    """

    def make_proximal_composition(step_size):
        """Create proximal for the dual with a given step_size

        Parameters
        ----------
        step_size : positive `float`
            Step size parameter

        Returns
        -------
        proximal : `Operator`
            The proximal operator of ``prox[step_size * f * L](x)``
        """
        Id = IdentityOperator(operator.domain)
        Ir = IdentityOperator(operator.range)
        prox_muf = proximal(step_size)
        return Id + (1.0 / mu) * operator.adjoint((prox_muf - Ir) * operator)

    return make_proximal_composition


def proximal_zero(space):
    """Function to create the proximal operator of the zero functional.

    Function to initialize the proximal operator of the zero functional
    defined on ``space``. The proximal operator of this functional is the
    identity operator

        prox[tau * G](x) = x  where G=0

    It is independent of tau.

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp` spaces
        Domain of the functional G=0

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized
    """

    def make_identity(tau):
        """Return an instance of the proximal operator.

        Parameters
        ----------
        tau : positive `float`
            Unused step size parameter. Introduced to provide a unified
            interface

        Returns
        -------
        id : `IdentityOperator`
            The proximal operator instance of G = 0 which is the
            identity operator
        """

        return IdentityOperator(space)

    return make_identity


def proximal_clamp(space, lower=None, upper=None):
    """Function to create the proximal operator of G(x) = ind(a <= x <= b).

    Function for the proximal operator of the functional G(x) = ind(a<=x<=b)
    to be initialized.

    If P is the set of elements with a < x < b, the indicator function of
    which is defined as

        ind(a < x < b) = {0 if x in P, infinity if x is not in P}

    with x being an element in ``space``.

    The proximal operator of ``tau * G^*`` where ``tau`` is a step size
    is the point-wise non-negativity thresholding of x

                              a if x < a,
         prox[tau * G](x) = { x if a <= x <= b
                              b if x > b

    It is independent of tau and invariant under a positive rescaling of G
    which leaves the indicator function as it stands.

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp`
        Domain of the functional G(x)
    lower : ``space.field`` element or ``space`` element-like, optional
        The lower bound. Default: `None`, interpreted as -infinity
    upper : ``space.field`` element or ``space`` element-like, optional
        The upper bound. Default: `None`, interpreted as infinity

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    See Also
    --------
    proximal_nonnegativity : Special case with ``lower=0, upper=infty``
    """

    # Convert element-likes if needed, also does some space checking
    if lower is not None and lower not in space and lower not in space.field:
        lower = space.element(lower)
    if upper is not None and upper not in space and upper not in space.field:
        upper = space.element(upper)

    if lower in space.field and upper in space.field:
        if lower > upper:
            raise ValueError('Invalid values, `lower` ({}) > `upper` ({}).'
                             ''.format(lower, upper))
        assert lower < upper

    class _ProxOpClamp(Operator):

        """The proximal operator."""

        def __init__(self, tau):
            """Initialize the proximal operator.

            Parameters
            ----------
            tau : positive `float`
                Step size parameter. Not used.
            """
            self.tau = float(tau)
            super().__init__(domain=space, range=space, linear=False)

        def _call(self, x, out):
            """Apply the operator to ``x`` and store the result in ``out``."""

            # Point-wise non-negativity thresholding: x if x > 0, else 0
            if lower is not None and upper is None:
                x.ufunc.maximum(lower, out=out)
            elif lower is None and upper is not None:
                x.ufunc.minimum(upper, out=out)
            elif lower is not None and upper is not None:
                x.ufunc.maximum(lower, out=out)
                out.ufunc.minimum(upper, out=out)
            else:
                out.assign(x)

    return _ProxOpClamp


def proximal_nonnegativity(space):
    """Function to create the proximal operator of G(x) = ind(x >= 0).

    Function for the proximal operator of the functional G(x)=ind(x >= 0) to be
    initialized.

    If P is the set of non-negative elements, the indicator function of
    which is defined as

        ind(x >= 0) = {0 if x in P, infinity if x is not in P}

    with x being an element in ``space``.

    The proximal operator of ``tau * F^*`` where ``tau`` is a step size
    is the point-wise non-negativity thresholding of x

         prox[tau * G](x) = {x if x >= 0, 0 if < 0}

    It is independent of tau and invariant under a positive rescaling of G
    which leaves the indicator function as it stands.

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp`
        Domain of the functional G(x)

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    See Also
    --------
    proximal_clamp : The underlying implementation
    """

    return proximal_clamp(space, lower=0)


def proximal_convexconjugate_l2(space, lam=1, g=None):
    """Proximal operator factory of the convex conjugate of the l2-norm.

    Function for the proximal operator of the convex conjugate of the
    functional F where F is the l2-norm

        F(x) =  lam ||x - g||_2

    with x and g elements in ``space``, scaling factor lam, and given data g.

    The convex conjugate, F^*, of F is given by

        F^*(y) = {0 if ||x-g|| < lam, infty else}

    The proximal operator of ``sigma * F^*`` where ``sigma`` is a step size
    is given by

        prox[sigma * F^*](y) = (y - sigma g) / (1 + sigma/lam)

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp`
        Domain of F(x)
    g : `DiscreteLpVector`
        An element in ``space``
    lam : positive `float`
        Scaling factor or regularization parameter

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    Notes
    -----
    Most problems are forumlated for the norm squared, in that case use that
    proximal operator instead.

    See Also
    --------
    proximal_l2 : proximal without convex conjugate
    proximal_convexconjugate_l2 : proximal without square
    """

    prox_l2 = proximal_l2(space, lam=lam, g=g)
    return proximal_convexconjugate(prox_l2)


def proximal_l2(space, lam=1, g=None):
    """Proximal operator factory of the l2-norm.

    Function for the proximal operator of the convex conjugate of the
    functional F where F is the l2-norm

        F(x) =  lam ||x - g||_2

    The proximal operator of ``sigma * F`` where ``sigma`` is a step size
    is given by

        prox[sigma * F^*](y) = (y - sigma g) / (1 + sigma/lam)

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp`
        Domain of F(x)
    g : `DiscreteLpVector`
        An element in ``space``
    lam : positive `float`
        Scaling factor or regularization parameter

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    Notes
    -----
    Most problems are forumlated for the norm squared, in that case use that
    proximal operator instead.

    See Also
    --------
    proximal_l2_squared : proximal for norm squared
    proximal_convexconjugate_l2 : proximal for convex conjugate
    """

    lam = float(lam)

    if g is not None and g not in space:
        raise TypeError('{!r} is not an element of {!r}'.format(g, space))

    class _ProximalL2(Operator):

        """The proximal operator."""

        def __init__(self, sigma):
            """Initialize the proximal operator.

            Parameters
            ----------
            sigma : positive `float`
                Step size parameter
            """
            self.sigma = float(sigma)
            super().__init__(domain=space, range=space, linear=False)

        def _call(self, x, out):
            """Apply the operator to ``x`` and stores the result in ``out``."""

            step = self.sigma * lam

            if g is None:
                x_norm = x.norm()

                if x_norm >= step:
                    out.lincomb(1.0 - step / x_norm, x)
                else:
                    out.set_zero()

            else:
                diff_norm = (x - g).norm()

                if diff_norm >= step:
                    out.lincomb(1.0 - step / diff_norm, x,
                                step / diff_norm, g)
                else:
                    out.assign(g)

    return _ProximalL2


def proximal_convexconjugate_l2_squared(space, lam=1, g=None):
    """Proximal operator factory of the convex conjugate of the l2-norm square.

    Function for the proximal operator of the convex conjugate of the
    functional F where F is the l2-norm

        F(x) =  lam 1/2 ||x - g||_2^2

    with x and g elements in ``space``, scaling factor lam, and given data g.

    The convex conjugate, F^*, of F is given by

        F^*(y) = 1/lam (1/2 ||y/lam||_2^2 + <y/lam,g>)

    The proximal operator of ``sigma * F^*`` where ``sigma`` is a step size
    is given by

        prox[sigma * F^*](y) = (y - sigma g) / (1 + sigma/lam)

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp`
        Domain of F(x)
    g : `DiscreteLpVector`
        An element in ``space``
    lam : positive `float`
        Scaling factor or regularization parameter

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    See Also
    --------
    proximal_convexconjugate_l2 : proximal without square
    proximal_l2_squared : proximal without convex conjugate
    """
    lam = float(lam)

    if g is not None and g not in space:
        raise TypeError('{!r} is not an element of {!r}'.format(g, space))

    class _ProximalConvConjL2Squared(Operator):

        """The proximal operator."""

        def __init__(self, sigma):
            """Initialize the proximal operator.

            Parameters
            ----------
            sigma : positive `float`
                Step size parameter
            """
            self.sigma = float(sigma)
            super().__init__(domain=space, range=space, linear=g is None)

        def _call(self, x, out):
            """Apply the operator to ``x`` and stores the result in
            ``out``"""

            # (x - sig*g) / (1 + sig/lam)

            sig = self.sigma
            if g is None:
                out.lincomb(1 / (1 + sig / lam), x)
            else:
                out.lincomb(1 / (1 + sig / lam), x, -sig / (1 + sig / lam), g)

    return _ProximalConvConjL2Squared


def proximal_l2_squared(space, lam=1, g=None):
    """Proximal operator factory of the l2-norm square.

    Function for the proximal operator of the convex conjugate of the
    functional F where F is the l2-norm

        F(x) =  lam 1/2 ||x - g||_2^2

    with x and g elements in ``space``, scaling factor lam, and given data g.

    The proximal operator of ``sigma * F^*`` where ``sigma`` is a step size
    is given by

        prox[sigma * F^*](y) = (y - sigma g) / (1 + sigma/lam)

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp`
        Domain of F(x)
    g : `DiscreteLpVector`
        An element in ``space``
    lam : positive `float`
        Scaling factor or regularization parameter

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    See Also
    --------
    proximal_l2 : proximal without square
    proximal_convexconjugate_l2_squared : proximal for convex conjugate
    """

    # TODO: optimize
    prox_cc_l2_squared = proximal_convexconjugate_l2_squared(space,
                                                             lam=lam, g=g)
    return proximal_convexconjugate(prox_cc_l2_squared)


def proximal_convexconjugate_l1(space, lam=1, g=None):
    """Proximal operator factory of the convex conjugate of the l1-semi-norm.

    Function for the proximal operator of the convex conjugate of the
    functional F where F is an l1-semi-norm

        F(x) = lam || ||x-g||_p ||_1

    with x and g elements in ``space``, scaling factor lam, and point-wise
    magnitude ||x||_p of x. If x is vector-valued, ||x||_p is the point-wise
    l2-norm across the vector components.

    The convex conjugate, F^*, of F is given by the indicator function of
    the set box(lam)

        F^*(y) = lam ind_{box(lam)}(||y / lam||_p + <y / lam, g>)

    where box(lam) is a hypercube centered at the origin with width 2 lam.

    The proximal operator of ``sigma * F^*`` where ``sigma`` is a step size
    is given by

        prox[sigma * F^*](y) =
            lam (y - sigma g) / (max(lam 1_{||y||_p}, ||y - sigma g||_p)

    where max(.,.) thresholds the lower bound of ||y||_p point-wise and
    1_{||y||_p} is a vector in the space of ||y||_p with all components set
    to 1.

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp` spaces
        Domain of the functional F
    g : `DiscreteLpVector`
        An element in ``space``
    lam : positive `float`
        Scaling factor or regularization parameter

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    See Also
    --------
    proximal_l1 : proximal without convex conjugate conjugate
    """
    lam = float(lam)

    if g is not None and g not in space:
        raise TypeError('{!r} is not an element of {!r}'.format(g, space))

    class _ProximalConvConjL1(Operator):

        """The proximal operator."""

        def __init__(self, sigma):
            """Initialize the proximal operator.

            Parameters
            ----------
            sigma : positive `float`
                Stepsize parameter
            """
            # sigma is not used
            self.sigma = float(sigma)
            super().__init__(domain=space, range=space, linear=False)

        def _call(self, x, out):
            """Apply the operator to ``x`` and stores the result in ``out``."""

            # lam * (x - sigma * g) / max(lam, |x - sigma * g|)

            if g is not None:
                diff = x - self.sigma * g
            else:
                diff = x

            if isinstance(x.space, ProductSpace):
                # Calculate |x| = pointwise 2-norm of x

                tmp = diff[0] ** 2
                sq_tmp = x[0].space.element()
                for x_i in diff[1:]:
                    sq_tmp.multiply(x_i, x_i)
                    tmp += sq_tmp
                tmp.ufunc.sqrt(out=tmp)

                # Pointwise maximum of |x| and lambda
                tmp.ufunc.maximum(lam, out=tmp)

                # Global scaling
                tmp /= (lam)

                # Pointwise division
                for out_i, x_i in zip(out, diff):
                    out_i.divide(x_i, tmp)

            else:
                # Calculate |x| = pointwise 2-norm of x
                diff.ufunc.absolute(out=out)

                # Pointwise maximum of |x| and lambda
                out.ufunc.maximum(lam, out=out)

                # Global scaling
                out /= (lam)

                # Pointwise division
                out.divide(diff, out)

    return _ProximalConvConjL1


def proximal_l1(space, lam=1, g=None):
    """Proximal operator factory of the l1-semi-norm.

    Function for the proximal operator of the functional F where F is an
    l1-semi-norm

        F(x) = lam || ||x-g||_p ||_1

    with x and g elements in ``space``, scaling factor lam, and point-wise
    magnitude ||x||_p of x. If x is vector-valued, ||x||_p is the point-wise
    l2-norm across the vector components.

    The proximal operator of F is

                              y - lam   if y > lam,
         prox[lam * F](y) = { 0         if -lam <= y <= lam
                              y + lam   if y < -lam


    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp` spaces
        Domain of the functional F
    g : `DiscreteLpVector`
        An element in ``space``
    lam : positive `float`
        Scaling factor or regularization parameter

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    See Also
    --------
    proximal_convexconjugate_l1 : proximal for convex conjugate
    """

    # TODO: optimize
    prox_cc_l1 = proximal_convexconjugate_l1(space, lam=lam, g=g)
    return proximal_convexconjugate(prox_cc_l1)


def proximal_convexconjugate_kl(space, lam=1, g=None):
    """Proximal operator factory of the convex conjugate of the KL divergence.

    Function returning the proximal operator of the convex conjugate of the
    functional F where F is the entropy-type Kullback-Leibler (KL) divergence

        F(x) = sum_i (x - g + g ln(g) - g ln(pos(x)))_i + ind_P(x)

    with x and g in X and g non-negative. The indicator function ind_P(x)
    for the positive elements of x is used to restrict the domain of F such
    that F is defined over whole X. The non-negativity thresholding pos is
    used to define F in the real numbers.

    The proximal operator of the convex conjugate, F^*, of F is

        F^*(p) = sum_i (-g ln(pos(1_X - p))_i + ind_P(1_X - p)

    where p is the variable dual to x, and 1_X is a vector in the space X with
    all components set to 1.

    The proximal operator of the convex conjugate of F is

        prox[sigma * F^*](x) =
            1/2 (lam + x - sqrt((x - lam)^2 + 4 lam sigma g)

    with the step size parameter sigma and lam_X is a vector in the space X
    with all components set to lam.

    Parameters
    ----------
    space : `DiscreteLp` or `ProductSpace` of `DiscreteLp` spaces
        The space X which is the domain of the functional F
    g : `DiscreteLpVector`
        An element in ``space``
    lam : positive `float`
        Scaling factor

    Returns
    -------
    prox : `callable`
        Factory for the proximal operator to be initialized

    Notes
    -----
    KL based objectives are common in MLEM optimization problems and are often
    used when data noise governed by a multivariate Poisson probability
    distribution is significant.

    The intermediate image estimates can have negative values even though
    the converged solution will be non-negative. Non-negative intermediate
    image estimates can be enforced by adding an indicator function ind_P
    the primal objective.
    """
    lam = float(lam)

    if g is not None and g not in space:
        raise TypeError('{} is not an element of {}'.format(g, space))

    class _ProximalConvConjKL(Operator):

        """The proximal operator."""

        def __init__(self, sigma):
            """Initialize the proximal operator.

            Parameters
            ----------
            sigma : positive `float`
            """
            self.sigma = float(sigma)
            super().__init__(domain=space, range=space, linear=False)

        def _call(self, x, out):
            """Apply the operator to ``x`` and stores the result in ``out``."""

            # 1 / 2 (lam_X + x - sqrt((x - lam_X) ^ 2 + 4; lam sigma g)

            # out = x - lam_X
            out.assign(x)
            out -= lam

            # (out)^2
            out.ufunc.square(out=out)

            # out = out + 4 lam sigma g
            if g is not None:
                out.lincomb(1, out, 4.0 * lam * self.sigma, g)

            # out = sqrt(out)
            out.ufunc.sqrt(out=out)

            # out = x - out
            out.lincomb(1, x, -1, out)

            # out = lam_X + out
            out.lincomb(lam, space.one(), 1, out)

            # out = 1/2 * out
            out /= 2

    return _ProximalConvConjKL


def proximal_variable_lp(space, exponent, lam=1.0, g=None):
    """Return the proximal operator of the variable Lebesgue modular.

    Parameters
    ----------
    space : `DiscreteLp` or power space of such
        Space on which the proximal operator acts
    exponent : ``space`` element-like or `float`
        Variable (or constant) exponent used in the modular
    lam : positive `float`
        Scaling factor or regularization parameter
    g : ``space`` element-like
        An element in ``space``

    Notes
    -----
    The variable :math:`L^p` modular is the integral

        :math:`\\rho_p(f) := \int_{\Omega}
        \\lvert f(x) \\rvert^{p(x)}\, \mathrm{d}x`

    for :math:`\Omega \subset \mathbb{R}^d`, a function
    :math:`f:\Omega \\to \mathbb{R}^m` and an exponent mapping
    :math:`p:\Omega \\to [0, \infty)`.
    In this functional, an additional element :math:`g` can be
    specified, such that instead of :math:`f`, the difference
    :math:`f - g` is considered.
    """

    class VarLpModularProx(Operator):

        """Proximal operator of the variable Lebesgue modular."""

        def __init__(self, sigma):
            """Initialize a new instance.

            Parameters
            ----------
            sigma : positive `float`
                Scaling parameter in the proximal operator

            Notes
            -----
            The proximal operator of the variable :math:`L^p` modular
            can be defined point-wise almost everywhere as

                :math:`\mathrm{prox}_{\sigma \\rho_p}(f)(x)
                = \mathrm{arg}\min_{v \\in \mathbb{R}^m}
                \\left[\\lvert v \\rvert^{p(x)} + \\frac{1}{2\sigma}
                \\lvert v - f(x) \\rvert^2 \\right]`.

            For those :math:`x` where :math:`p(x) = 2`, we have

                :math:`\mathrm{prox}_{\sigma \\rho_p}(f)(x)
                = \\frac{f(x)}{1 + 2\sigma}`.

            In points with :math:`p(x) = 1`, it is

                :math:`\mathrm{prox}_{\sigma \\rho_p}(f)(x)
                = \max\\left\{1 - \\frac{\sigma}{\\lvert f(x) \\rvert},
                \, 0\\right\}\, f(x)`.

            Otherwise, the minimization problem is solved by a Newton
            method for p other than 1 or 2.
            """
            # TODO: check that we have a power space
            super().__init__(domain=space, range=space, linear=False)
            self.sigma = float(sigma)

            if isinstance(self.domain, ProductSpace):
                base_space = self.domain[0]
            else:
                base_space = self.domain

            if np.isscalar(exponent):
                self.exponent = float(exponent) * base_space.one()
            else:
                self.exponent = base_space.element(exponent)

            if g is not None:
                self.g = self.domain.element(g)
            else:
                self.g = None

        def _newton_iter(self, it, val, sigma, p, pm1, pm2, niter=5,
                         start_relax=0.5, tmp=None):
            """Helper method for the inner Newton iteration."""
            # Start value which guarantees convergence.
            np.divide(val, p, out=it)
            it /= pm2
            it /= -sigma
            # Avoid huge values by setting the start value to 1 if the
            # computed one is larger than 1. Since the power 1/(p-1) can be
            # as large as 20, we don't want large arguments. In those
            # points, a very small result is likely anyway.
            small_enough = (it <= 1)
            it[small_enough] **= 1.0 / pm1[small_enough]
            it[~small_enough] = 1.0
            it *= start_relax

            # The iteration itself
            for i in range(niter):
                # Denominator 1 + p*(p-1)*sigma*q**(p-2)
                tmp = np.power(it, pm2, out=tmp)
                tmp *= p
                tmp *= pm1
                tmp *= sigma
                tmp += 1.0

                # Nominator p*(p-2)*sigma*q**(p-1) + val
                np.power(it, pm1, out=it)
                it *= p
                it *= pm2
                it *= sigma
                it += val

                it /= tmp

        def _call(self, f, out, **kwargs):
            """Implement ``self(x, out, **kwargs)``.

            Parameters
            ----------
            f : domain element
                Element at which to evaluate the operator
            out : range element
                Element to which the result is written
            max_newton_iter : `int`, optional
                Maximum number of Newton iterations
            """
            if isinstance(self.domain, ProductSpace):
                self._call_pspace(f, out, **kwargs)
            else:
                self._call_scalar(f, out, **kwargs)

        def _call_scalar(self, f, out, **kwargs):
            """Implement ``self(x, out, **kwargs)`` for scalar domain."""
            if self.g is not None:
                f = f - self.g

            step = self.sigma * float(lam)

            exp_arr = self.exponent.asarray()
            out_arr = out.asarray()
            f_arr = f.asarray()
            f_nz = (f_arr != 0)

            # p = 2
            # This formula is used globally since it sets out to 0
            # where f is 0.
            out.lincomb(0, out, 1.0 / (1.0 + 2.0 * step), f)

            # p = 1 (taking also close to one for stability)
            cur_exp = (exp_arr <= 1.05)
            current = cur_exp & f_nz
            factor = np.maximum(1.0 - step / np.abs(f_arr[current]), 0.0)
            out_arr[current] = factor * f_arr[current]

            # Newton iteration for other p values. We consider only those
            # entries that correspond to f != 0.
            cur_exp = ~((exp_arr >= 1.95) | cur_exp)
            current = cur_exp & f_nz
            exp_p = exp_arr[current]
            exp_m1 = exp_p - 1
            exp_m2 = exp_p - 2
            it = out_arr[current]
            val = f_arr[current]
            tmp = np.empty_like(it)

            maxiter = int(kwargs.pop('max_newton_iter', 5))
            self._newton_iter(it, np.abs(val), step, exp_p, exp_m1, exp_m2,
                              niter=maxiter, tmp=tmp)

            out_arr[current] = it
            out_arr[current] *= np.sign(val)
            out[:] = out_arr

            if self.g is not None:
                out += self.g

        def _call_pspace(self, f, out, **kwargs):
            """Implement ``self(x, out, **kwargs)`` for vectorial domain."""
            if self.g is not None:
                f = f - self.g

            step = self.sigma * float(lam)

            exp_arr = self.exponent.asarray()
            f_nz = [(fi.asarray() != 0) for fi in f]
            pw_norm = PointwiseNorm(self.domain)
            f_norm = pw_norm(f)

            # p = 2
            # This formula is used globally since it sets out to 0
            # where f is 0.
            out.lincomb(0, out, 1.0 / (1.0 + 2.0 * step), f)

            # p = 1 (taking also close to one for stability)
            cur_exp = (exp_arr <= 1.05)
            for fi, fi_nz, oi in zip(f, f_nz, out):
                fi_arr = fi.asarray()
                oi_arr = oi.asarray()
                current = cur_exp & fi_nz
                factor = np.maximum(1.0 - step / np.abs(fi_arr[current]), 0.0)
                oi_arr[current] = factor * fi_arr[current]
                oi[:] = oi_arr

            # Newton iteration for other p values
            maxiter = int(kwargs.pop('max_newton_iter', 5))
            cur_exp = ~((exp_arr >= 1.95) | cur_exp)
            for fi, fi_nz, oi in zip(f, f_nz, out):
                fi_arr = fi.asarray()
                oi_arr = oi.asarray()
                current = cur_exp & fi_nz
                exp_p = exp_arr[current]
                exp_m1 = exp_p - 1
                exp_m2 = exp_p - 2

                it = oi_arr[current]
                val = fi_arr[current]
                tmp = np.empty_like(it)
                self._newton_iter(it, np.abs(val), step, exp_p, exp_m1, exp_m2,
                                  niter=maxiter, tmp=tmp)

                oi_arr[current] = it
                oi_arr[current] /= f_norm.asarray()[current]
                oi_arr[current] *= val

                oi[:] = oi_arr

            if self.g is not None:
                out += self.g

    return VarLpModularProx


if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
