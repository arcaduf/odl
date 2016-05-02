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

import numpy as np

from odl.discr.discr_ops import Gradient, Divergence
from odl.operator.default_ops import ResidualOperator
from odl.operator.operator import Operator, OpDomainError, OpRangeError
from odl.solvers.scalar.steplen import (
    StepLength, StepLengthFromIter, ConstantStepLength)
from odl.space.pspace import ProductSpace

def primal_dual_hybrid_gradient_tv(fwd_op, data, primal, dual, pstep, dstep,
                                   reg_param, proj_op, niter=1, callback=None):

    """Primal-dual method for TV regularisation.

    This is an implementation of the classical PDHG method by Zhu
    and Chan [ZC2008]_. It aims to minimise

        ``1/2 * ||T(f) - g||^2 + lambda * ||grad(f)||``

    where the first norm is a Hilbert space norm and the second one
    an arbitrary vector field norm. See Notes for mathematical details.

    Parameters
    ----------
    fwd_op : `Operator`
        Forward operator ``T`` of the reconstruction problem
    data : ``fwd_op.range`` `element-like`
        Data ``g`` from which to reconstruct the unknown ``f``
    primal : ``fwd_op.domain`` `element`
        Variable ``f`` used in the primal iteration. The final
        reconstruction will be stored in this object when the method
        is finished.
    dual : ``ProductSpace(primal.space, ndim)`` `element`
        An element of ``X^d``, where ``X`` is the reconstruction
        space and ``d`` is the number of dimensions. The space of
        this object is the same as the range of
        ``Gradient(primal.space)``. The dual iteration is performed
        in-place in this variable.
    pstep, dstep : {`float`, `iterator`, `StepLength`}
        Rule for the selection of the primal/dual step parameter
        ``theta/tau``.
        A float is interpreted as constant step, an iterator is queried
        in each iteration for a new value, and a `StepLength` object
        can be used for an adaptive rule.
    reg_param : nonnegative `float`
        Regularisation parameter ``lambda`` steering the weight of the
        total variation penaliser. Larger value means less trust in the
        data and more weight on the regulariser.
    proj_op : `Operator`
        Projection operator applied to the dual in each step. A proper
        implementation of this step is crucial since it is the only
        place where the norm of the gradient enters the iteration.
        The projection operator must map the ``dual.space`` to
        itself.
    niter : positive `int`, optional
        Number of iterations to compute
    callback : `callable`, optional
        Function or object to be evaluated on the primal in each step

    See also
    --------
    StepLength : rules for adaptive step length computation

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
    # --- Handle input parameters --- #

    if not isinstance(fwd_op, Operator):
        raise TypeError('{!r} is not an Operator instance.'.format(fwd_op))

    data = fwd_op.range.element(data)

    if primal not in fwd_op.domain:
        raise TypeError('primal {!r} is not an element of the forward '
                        'operator domain {!r}.'
                        ''.format(primal, fwd_op.domain))

    dual_space = ProductSpace(primal.space, primal.space.ndim)
    if dual not in dual_space:
        raise TypeError('dual {!r} is not an element of the power space '
                        'X^{} with X = primal space {!r}.'
                        ''.format(dual, primal.space.ndim, primal.space))

    # TODO: StepLength currently coexists with LineSearch, and the new
    # StepLengthFromIter is the only implementation.
    if np.isscalar(pstep):
        pstep = ConstantStepLength(pstep)
    elif isinstance(pstep, StepLength):
        pass
    else:
        pstep = StepLengthFromIter(iter(pstep))

    if np.isscalar(dstep):
        dstep = ConstantStepLength(dstep)
    elif isinstance(dstep, StepLength):
        pass
    else:
        dstep = StepLengthFromIter(iter(dstep))

    if reg_param < 0:
        raise ValueError('expected nonnegative regularization parameter, '
                         'got {}.'.format(reg_param))

    if not isinstance(proj_op, Operator):
        raise TypeError('{!r} is not an Operator instance.'.format(proj_op))
    if proj_op.domain != dual_space:
        raise OpDomainError('projection operator domain is not the dual '
                            'variable space {!r}.'.format(dual_space))
    if proj_op.range != dual_space:
        raise OpDomainError('projection operator range is not the dual '
                            'variable space {!r}.'.format(dual_space))

    # --- Create auxiliary objects --- #

    # TODO: generalise to arbitrary sparsifying operator instead of the
    # gradient?
    grad_op = Gradient(primal.space)
    # TODO: Why the hell is the divergence initialized with the primal space?
    # See https://github.com/odlgroup/odl/issues/375
    div_op = Divergence(primal.space)
    res_op = ResidualOperator(fwd_op, data)
    tmp_primal = primal.space.element()
    tmp_dual = dual.space.element()
    tmp_data = fwd_op.range.element()

    for _ in range(niter):
        # Dual step
        # TODO: define dstep interface, currently ignores input
        cur_dstep = dstep(dual)
        tmp_dual.lincomb(1, dual, -cur_dstep, grad_op(primal))
        proj_op(tmp_dual, out=dual)

        # Primal step
        cur_pstep = pstep(primal)
        res_op(primal, out=tmp_data)
        # Doing the evaluation of the update parts for the primal in two
        # steps to save one temporary
        # TODO: extension to non-linear operators?
        fwd_op.adjoint(tmp_data, out=tmp_primal)
        primal.lincomb(1, primal, -cur_pstep, tmp_primal)
        div_op(dual, out=tmp_primal)
        primal.lincomb(1, primal, -cur_pstep * reg_param, tmp_primal)
