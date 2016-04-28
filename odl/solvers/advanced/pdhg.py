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

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()


def primal_dual_hybrid_gradient_tv():

    """Primal-dual method for TV regularisation.

    This is an implementation of the classical PDHG method by Zhu
    and Chan [ZC2008]_. See Notes for mathematical details.

    Notes
    -----

    The PDHG method solves for given data :math:`g\in\mathcal{Y}`,
    a linear operator :math:`\mathcal{T}:\mathcal{X}\\to\mathcal{Y}`
    and a regularisation parameter :math:`\lambda > 0` the problem

        :math:`\min_{f \\in \mathcal{X}} \Big(\\frac{1}{2}
        \\lVert \mathcal{T}(f) - g\\rVert_{\mathcal{Y}}^2 +
        \\lambda \\lVert \mathrm{grad}f\\rVert \Big)`,

    where the first norm is the Hilbert space norm in the data space
    and the second a vectorial norm in the reconstruction space. We
    assume that the norm in the regularising term satisfies

        :math:`\\lVert f\\rVert =
        \\sup_{\\lVert h\\rVert_\\ast \\leq 1}
        \\langle f, h\\rangle_{\mathcal{X}}`,

    with the *associate norm* :math:`\\lVert \cdot \\rVert_\\ast`.
    Then the method can be written as the following iteration scheme in
    the primal variable :math:`p` and the dual variable :math:`d`:

        :math:`\mathbf{d}_{k+1} =
        \mathcal{P}_{\\lVert\cdot\\rVert_\\ast \\leq 1}(\mathbf{d}_k -
        \\tau\, \mathrm{grad} f)`,

        :math:`p_{k+1} =
        p_k - \\theta \\big(\mathcal{T}^*\\big(\mathcal{T}(p_k)
        - g\\big) + \\lambda\, \mathrm{div}\mathbf{d}_{k+1} \\big)`.

    Here, :math:`\mathcal{P}_{\\lVert\cdot\\rVert_\\ast \\leq 1}` denotes
    the projection onto the unit ball in the associate norm.
    """