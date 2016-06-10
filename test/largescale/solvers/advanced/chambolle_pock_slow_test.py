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

"""Tests for the chambolle pock method."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

# External
import numpy as np
import pytest

# Internal
import odl
from odl.util.testutils import example_element, almost_equal


operator_params = ['id', 'square']
operator_ids = [' operator = {} '.format(p) for p in operator_params]


@pytest.fixture(scope="module", ids=operator_ids, params=operator_params)
def operator(request):

    if request.param == 'id':
        space = odl.uniform_discr(0, 1, 5)
        return odl.IdentityOperator(space)

    elif request.param == 'square':
        space = odl.uniform_discr(0, 1, 5)
        return odl.PowerOperator(space, 2)


def test_chambolle_pock_inversion(operator):
    # Starting point for explicit computation
    discr_vec_0 = example_element(operator.domain)

    # Copy to be overwritten by the algorithm
    discr_vec = discr_vec_0.copy()

    # Proximal operator using the same factory function for F^* and G
    prox_dual = odl.solvers.proximal_cconj_l2_squared(operator.domain)
    prox_primal = odl.solvers.proximal_zero(operator.range)

    opnorm = odl.power_method_opnorm(operator, 10, xstart=discr_vec)
    tau = sigma = 0.5 / opnorm ** 2

    # Run the algorithm
    odl.solvers.chambolle_pock_solver(operator, discr_vec,
                                      tau=tau, sigma=sigma,
                                      proximal_primal=prox_primal,
                                      proximal_dual=prox_dual,
                                      niter=1000)

    print(discr_vec)
    print(discr_vec_0)
    print(tau, opnorm)

    assert almost_equal(discr_vec.norm(), 0.0, places=2)


if __name__ == '__main__':
    pytest.main(str(__file__.replace('\\', '/') + ' -v --largescale'))
