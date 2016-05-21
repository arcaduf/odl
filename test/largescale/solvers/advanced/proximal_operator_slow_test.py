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

"""Tests for the factory functions to create proximal operators."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

# External
import numpy as np
import pytest

# Internal
import odl
from odl.solvers.advanced.proximal_operators import (
    combine_proximals, proximal_zero, proximal_nonnegativity,
    proximal_l1, proximal_convexconjugate_l1,
    proximal_l2, proximal_convexconjugate_l2,
    proximal_l2_squared, proximal_convexconjugate_l2_squared,
    proximal_convexconjugate_kl, proximal_variable_lp)
from odl.util.testutils import example_element


pytestmark = odl.util.skip_if_no_largescale


class VariableLpModular(odl.Operator):

    """Functional for evaluating the variable Lp modular.

    The variable Lp modular is defined as

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
            Discretized function space on which the modular is defined
        var_exp : scalar-valued ``space`` `element-like`
            The variable exponent ``p(x)``
        """
        super().__init__(space, space.field, linear=False)
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


step_params = [0.1, 1.0, 10.0]
step_ids = [' stepsize = {} '.format(p) for p in step_params]


@pytest.fixture(scope="module", ids=step_ids, params=step_params)
def stepsize(request):
    return request.param


offset_params = [False, True]
offset_ids = [' offset = {} '.format(str(p).ljust(5)) for p in offset_params]


@pytest.fixture(scope="module", ids=offset_ids, params=offset_params)
def offset(request):
    return request.param


prox_params = ['l1 ', 'l1_dual',
               'l2', 'l2_dual',
               'l2^2', 'l2^2_dual',
               'var_lp_mod']
prox_ids = [' f = {}'.format(p.ljust(10)) for p in prox_params]


@pytest.fixture(scope="module", ids=prox_ids, params=prox_params)
def proximal_and_function(request):
    """Return a proximal factory and the corresponding function."""
    name = request.param.strip()

    space = odl.uniform_discr(0, 1, 2)

    if name == 'l1':
        def l1_norm(x):
            return np.abs(x).inner(x.space.one())

        prox = proximal_l1(space)

        return prox, l1_norm

    if name == 'l1_dual':
        def l1_cc_norm(x):
            sup_norm = np.max(np.abs(x))
            if sup_norm <= 1.0:
                return 0.0
            else:
                return np.Infinity

        prox = proximal_convexconjugate_l1(space)

        return prox, l1_cc_norm

    elif name == 'l2':
        def l2_norm(x):
            return x.norm()

        prox = proximal_l2(space)

        return prox, l2_norm

    elif name == 'l2_dual':
        def l2_norm_dual(x):
            if x.norm() <= 1.00001:  # numerical margin
                return 0.0
            else:
                return np.Infinity

        prox = proximal_convexconjugate_l2(space)

        return prox, l2_norm_dual

    elif name == 'l2^2':
        def l2_norm_squared(x):
            return 0.5 * x.norm() ** 2

        prox = proximal_l2_squared(space)

        return prox, l2_norm_squared

    elif name == 'l2^2_dual':
        def l2_norm_squared(x):
            return 0.5 * x.norm() ** 2

        prox = proximal_convexconjugate_l2_squared(space)

        return prox, l2_norm_squared

    elif name == 'var_lp_mod':
        exponent = np.random.uniform(1, 2, space.shape)
        exponent[exponent <= 1.01] = 1.0
        exponent[exponent >= 1.99] = 2.0
        var_lp_mod = VariableLpModular(space, exponent)
        prox = proximal_variable_lp(space, exponent)
        return prox, var_lp_mod

    else:
        assert False


def proximal_objective(function, step_size, x, y):
    """Calculate the objective function of the proximal optimization problem"""
    return function(y) + 1.0 / (2.0 * step_size) * (x - y).norm() ** 2


def test_proximal_defintion(proximal_and_function, stepsize):
    """Test the defintion of the proximal:

        prox[lam * f](x) = argmin_y {f(y) + 1/(2 lam) ||x-y||}

    Hence we expect for all x in the domain of the proximal

        x* = prox[lam * f](x)
        f(x*) + 1/(2 lam) ||x*-y|| < f(y) + 1/(2 lam) ||x-y||
    """

    proximal_factory, function = proximal_and_function

    proximal = proximal_factory(stepsize)

    assert proximal.domain == proximal.range

    x = example_element(proximal.domain) * 10
    f_x = proximal_objective(function, stepsize, x, x)
    prox_x = proximal(x)
    f_prox_x = proximal_objective(function, stepsize, x, prox_x)

    print(prox_x, x)
    print(prox_x.norm(), x.norm())
    print(f_prox_x, f_x)
    assert f_prox_x <= f_x

    for i in range(100):
        y = example_element(proximal.domain)
        f_y = proximal_objective(function, stepsize, x, y)

        print(prox_x, y)
        print(f_prox_x, f_y)
        assert f_prox_x <= f_y

if __name__ == '__main__':
    pytest.main(str(__file__.replace('\\', '/') + ' -v --largescale'))
