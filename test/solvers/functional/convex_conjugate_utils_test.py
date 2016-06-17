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

"""Tests for the utility functions for utilities/algebra for convex conjugate
functionals."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

# External
import numpy as np
import pytest

# Internal
import odl
import odl.solvers as odls

from odl.util.testutils import all_almost_equal, example_element

# Places for the accepted error when comparing results
PLACES = 8


# TODO: Remove this one when the "Functional" branch is in
class L2NormSquare(odls.Functional):
    def __init__(self, domain):
        super().__init__(domain=domain, linear=False, convex=True,
                         concave=False, smooth=True, grad_lipschitz=2)

    def _call(self, x):
        return np.abs(x).inner(np.abs(x))

    @property
    def gradient(self):
        functional = self

        class L2SquareGradient(odl.Operator):
            def __init__(self):
                super().__init__(functional.domain, functional.domain,
                                 linear=False)

            def _call(self, x):
                return 2.0*x

        return L2SquareGradient()

    @property
    def conjugate_functional(self):
        functional = self

        class L2SquareConjugateFunctional(odls.Functional):
            def __init__(self):
                super().__init__(functional.domain, linear=False)

            def _call(self, x):
                return x.norm()**2 * 1/4

            @property
            def gradient(self):
                functional = self

                class L2CCSquareGradient(odl.Operator):
                    def __init__(self):
                        super().__init__(functional.domain, functional.domain,
                                         linear=False)

                    def _call(self, x):
                        return x*(1.0/2.0)

                return L2CCSquareGradient()

        return L2SquareConjugateFunctional()


def test_convex_conjugate_translation():
    """Test for the convex conjugate of a translation: (f(. - y))^*"""

    # Image space
    space = odl.uniform_discr(0, 1, 10)

    # Element in the image space where the proximal operator is evaluated
    translation = example_element(space)

    # Creating the functional ||x||_2^2
    test_functional = L2NormSquare(space)
    cc_test_functional = test_functional.conjugate_functional

    # Testing that translation belonging to the wrong space gives TypeError
    wrong_space = odl.uniform_discr(1, 2, 10)
    wrong_translation = example_element(wrong_space)
    with pytest.raises(TypeError):
        odls.convex_conjugate_translation(cc_test_functional,
                                          wrong_translation)

    # Create translated convex conjugate functional
    cc_translated = odls.convex_conjugate_translation(cc_test_functional,
                                                      translation)

    # Create an element in the space, in which to evaluate
    x = example_element(space)

    # Test for evaluation of the functional
    # Explicit computation: 1/4 * ||x||^2 + <x,translation>
    expected_result = 1.0/4.0 * x.norm()**2 + x.inner(translation)
    assert all_almost_equal(cc_translated(x), expected_result, places=PLACES)

    # Test for the gradient
    # Explicit computation: x/2 + translation
    expected_result = 1.0/2.0 * x + translation
    cc_translated_gradient = cc_translated.gradient
    assert all_almost_equal(cc_translated_gradient(x), expected_result,
                            places=PLACES)

    # Test for derivative in direction p
    p = example_element(space)

    # Explicit computation in point x, in direction p: <x/2 + translation, p>
    expected_result = p.inner(1.0/2.0 * x + translation)
    assert all_almost_equal(cc_translated.derivative(x)(p), expected_result,
                            places=PLACES)

if __name__ == '__main__':
    pytest.main(str(__file__.replace('\\', '/') + ' -v'))
