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

"""Example showing how to use vectorization of FunctionSpaceVector's."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

import numpy as np
import odl


class L2DataMatchingFunctional(odl.solvers.Functional):

    """Basic data matching functional using the L2 norm.

    This functional computes::

        1/2 * ||f - g||_2^2

    for a given element ``g``.
    """

    def __init__(self, space, data):
        """Initialize a new instance.

        Parameters
        ----------
        space : `DiscreteLp` with exponent 2.0
            Space where the data is matched
        data : `DiscreteLp` element-like
            Data which is to be matched
        """
        if not (isinstance(space, odl.DiscreteLp) and space.exponent == 2.0):
            raise ValueError('not an L2 space.')
        super().__init__(space, linear=False)
        self.data = self.domain.element(data)

    def _call(self, x):
        """Return ``self(x)``."""
        return self.domain.dist(x, self.data)

    def gradient(self, x):
        """Return the gradient in the point ``x``."""
        return x - self.data

    def derivative(self, x):
        """Return the derivative in ``x``."""
        return self.gradient(x).T


class VariableLpModulus(odl.solvers.Functional):

    """Functional for evaluating the variable Lp modulus.

    The variable Lp modulus is defined as

        ``S(f) = integral( |f(x)|^(p(x)) dx ) + ||f_inf||_inf``

    where ``f_inf`` is the restriction of ``f`` to the set where
    ``p = inf`` and the integral is taken over the set where ``p``
    is finite.

    It maps a real-valued function on Omega to the real numbers.
    The exponent function is expected to be fulfill
    ``1 <= p(x) <= inf`` everywhere.
    """

    def __init__(self, space, var_exp):
        """Initialize a new instance.

        Parameters
        ----------
        space : `DiscreteLp`
            Discretized function space on which the modulus is defined
        var_exp : scalar-valued ``space`` `element-like`
            The variable exponent ``p(x)``
        """
        super().__init__(space, linear=False)
        self.var_exp = self.domain.element(var_exp)
        self._var_exp_flat = self.var_exp.asarray().ravel()
        self._fin_idcs = np.where(np.isfinite(self._var_exp_flat))
        self._inf_idcs = np.where(np.isinf(self._var_exp_flat))
        self._all_finite = (self._inf_idcs[0].size == 0)

    def _call(self, f):
        """Return ``self(f)``."""
        if self._all_finite:
            tmp = np.power(np.abs(f), self.var_exp)
            return self.domain.inner(tmp, self.domain.one())

        else:
            # Integral part
            tmp = self.domain.zero()
            absf = np.abs(f)
            tmp[self._fin_idcs] = np.power(absf[self._fin_idcs],
                                           self.var_exp[self._fin_idcs])
            fin_term = self.domain.inner(tmp, self.domain.one())

            # Inf part
            tmp[:] = 0
            tmp[self._inf_idcs] = absf[self._inf_idcs]
            inf_term = np.amax(tmp)

            return fin_term + inf_term

    def gradient(self, f, out=None):
        """Evaluate the gradient in ``x``.

        The gradient is given as

            ``grad S(f) = p * |f|^(p - 2) * f``

        in the sense of point-wise operations.
        """
        if not self._all_finite:
            raise NotImplementedError('gradient not well-defined for infinite '
                                      'exponent.')
        if out is None:
            out = self.domain.element()

        f.ufunc.absolute(out=out)
        out.ufunc.power(self.var_exp - 2, out=out)
        out.asarray()[np.isnan(out.asarray())] = 0  # Handle NaN
        out *= f
        out *= self.var_exp
        return out


class VariableLpNorm(odl.solvers.Functional):

    """The p-norm with spatially varying exponent ``p``.

    The variable Lp norm is defined as

        ``||f||_p = inf{s > 0 | rho_p(f / s) <= 1}``

    where ``rho_p`` is the variable Lp modulus. Starting from the
    initial guess ``s = rho_p(f)``, a bisection method is used to
    determine the optimal ``s``.
    """

    def __init__(self, space, var_exp):
        """Initialize a new instance.

        Parameters
        ----------
        space : `DiscreteLp`
            Discretized function space on which the modulus is defined
        var_exp : scalar-valued ``space`` `element-like`
            The variable exponent ``p(x)``
        """
        super().__init__(space, linear=False)
        self.var_exp = self.domain.element(var_exp)
        self._min_exp = np.min(self.var_exp)
        self.modulus = VariableLpModulus(space, var_exp)

    def _call(self, f, **kwargs):
        """Return ``self(f)``.

        Parameters
        ----------
        f : `DiscreteLpVector`
            Function whose norm to calculate
        atol : positive `float`, optional
            Stop the iteration in the norm computation when
            ``|rho_p(f / s) - 1| <= atol``.
            Default: 0.01
        maxiter : positive `int`, optional
            Iterate at most this many times. Default: 10
        """
        atol = kwargs.pop('atol', 0.01)
        maxiter = kwargs.pop('maxiter', 10)

        s = self.modulus(f)
        if s == 0:
            return 0.0

        m = self.modulus(f / s)
        if abs(m - 1) <= atol:
            return s
        elif m < 1:
            fac = 0.5
        else:
            fac = 2.0

        # Find a starting point for the s iteration
        m_old = m
        s_old = s
        it = 0
        while True:
            s *= fac
            m = self.modulus(f / s)
            it += 1
            if np.sign(m - 1) != np.sign(m_old - 1):
                break
            else:
                m_old = m
                s_old = s

        # Iterate until tolerance or maximum number of iterations reached
        s_low, s_up = min(s, s_old), max(s, s_old)
        for _ in range(maxiter - it + 1):
            s_test = (s_low + s_up) / 2  # TODO: use golden ratio
            m_test = self.modulus(f / s_test)
            if abs(m_test - 1) <= atol:
                return s_test
            elif m_test < 1:
                s_up = s_test
            else:
                s_low = s_test

        return (s_low + s_up) / 2


class VariableLpUnitBallProjector(odl.Operator):

    """Projector onto the unit ball ``{f: ||f||_p <= 1}``.

    Currently, we implement the simplified version

        ``P(f) = f / ||f||_p``.
    """

    def __init__(self, norm_func):
        """Initialize a new instance.

        Parameters
        ----------
        norm_func : `VariableLpNorm`
            Functional to evaluate the norm
        """
        self.norm = norm_func
        super().__init__(self.norm.domain, self.norm.domain, linear=False)

    def _call(self, x):
        """Return ``self(x)``."""
        return x / self.norm(x)


class GradientOperator(odl.Operator):

    """Wrap the gradient as operator to make the solvers happy."""

    def __init__(self, functional):
        assert isinstance(functional, odl.solvers.Functional)
        self._functional = functional
        super().__init__(self._functional.domain, self._functional.domain,
                         linear=self._functional.is_linear)

    def _call(self, x, out=None):
        return self._functional.gradient(x, out)


# ---- Example ---- #


def kernel(x):
    scaled_x = [xi / (np.sqrt(2) * 0.02) for xi in x]
    sq_sum = sum(xi ** 2 for xi in scaled_x)
    return np.exp(-sq_sum) / (2 * np.pi * 0.02 ** 2)


discr = odl.uniform_discr([-1, -1], [1, 1], (500, 500))
small_discr = odl.uniform_discr([-1, -1], [1, 1], (300, 300))
tmp = odl.util.phantom.submarine_phantom(small_discr, smooth=False)
tmp += odl.util.phantom.submarine_phantom(small_discr, smooth=True)
phantom = discr.zero()
phantom.asarray()[100:400, 100:400] = tmp
phantom.show('Phantom')

# Define the exponent: we use a convolution with a 3x3 constant kernel to
# broaden the boundaries (safety margin) without tapering off. After
# that, we set areas below a threshold to 2.0, above to 1.0.
# Possible issues:
# - Small jumps have smaller Laplacian -> low threshold needed, otherwise miss
# - Threshold is a parameter. How to coose?
# - Sensitivity to noise
# - Binary image - perhaps better to have a smooth function

exp_kernel = np.ones((3, 3))
exp_conv = odl.Convolution(discr, exp_kernel, impl='scipy_convolve',
                           scale=False)
lapl = odl.Laplacian(discr)
abs_lapl = np.abs(lapl(phantom))
abs_lapl /= np.amax(abs_lapl)
var_exponent = 2.0 - np.greater(exp_conv(abs_lapl), 0.3)
var_exponent.show('Exponent function')

# Set up the variable Lp TV functional - a composition
# lp_modulus o pointwise_norm o gradient
spatial_grad = odl.Gradient(discr)
pw_norm = odl.PointwiseNorm(spatial_grad.range)
var_lp_modulus = VariableLpModulus(discr, var_exp=var_exponent)
var_lp_tv_func = var_lp_modulus * pw_norm * spatial_grad

# Generate some data - a convolution of the phantom with a Gaussian kernel
conv = odl.Convolution(discr, kernel)
data = conv(phantom)
kernel_elem = discr.element(kernel)
kernel_elem.show('Convolution kernel')
data.show('Data - convolved phantom')


def simple_pdhg_norm(fwd_op, data, primal, dual, pstep, dstep, reg_param,
                     norm_op, maxiter=1):
    # Perform alternating steps in the primal and dual variables:
    #
    # * d_{k+1} = P(d_k - dstep * grad(p_k))
    # * p_{k+1} = p_k - pstep * (T^*(T(f) - g) + lambda * div(d_{k+1}))
    #
    # Here, P is a projector to the set {f: norm_op(f) <= 1}, which we
    # simplify to P(f) = f / norm_op(f), T is the (linear) forward operator,
    # and g is the data.
    # Note that the norm operator needs to be applied to the dual variables.
    grad_op = odl.Gradient(primal.space)
    div_op = odl.Divergence(primal.space)
    res_op = odl.ResidualOperator(fwd_op, data)
    tmp_primal_fwd = primal.space.element()
    tmp_primal_div = primal.space.element()
    tmp_data = fwd_op.range.element()

    for _ in range(maxiter):
        # Dual step
        dual.lincomb(1, dual, -dstep, grad_op(primal))
        dual /= norm_op(dual)

        # Primal step
        res_op(primal, out=tmp_data)
        fwd_op.adjoint(tmp_data, out=tmp_primal_fwd)
        div_op(dual, out=tmp_primal_div)
        primal.lincomb(1, primal, -pstep, tmp_primal_fwd)
        primal.lincomb(1, primal, -pstep * reg_param, tmp_primal_div)
        odl.util.numerics.apply_on_boundary(primal, lambda x: 0,
                                            out=primal.asarray())


# Vectorial variable exponent norm for the conjugate exponent
conj_max = 10.0  # values above this are considered infinite
small = np.where(var_exponent.asarray().ravel() <= 1 + 1 / (conj_max - 1.0))
conj_exp = discr.one() * np.inf
conj_exp[small] = 1.0 + 1.0 / (var_exponent[small] - 1)
norm_op = VariableLpNorm(discr, conj_exp) * pw_norm

# Step sizes
pstep = 0.2
dstep = 2.0

# Start values
primal = pstep * conv.adjoint(data)
odl.util.numerics.apply_on_boundary(primal, lambda x: 0, out=primal.asarray())
dual = spatial_grad.range.zero()
simple_pdhg_norm(conv, data, primal, dual, pstep, dstep, reg_param=0.01,
                 norm_op=norm_op, maxiter=2)
