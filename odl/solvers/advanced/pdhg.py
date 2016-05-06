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

from odl.operator.default_ops import ResidualOperator
from odl.operator.operator import Operator, OpDomainError, OpRangeError
from odl.solvers.scalar.steplen import (
    StepLength, StepLengthFromIter, ConstantStepLength)


__all__ = ('primal_dual_hybrid_gradient',)


def primal_dual_hybrid_gradient(fwd_op, data, primal, dual, spars_op,
                                proj_op, pstep, dstep, reg_param, niter=1,
                                **kwargs):
    """Primal-dual method for TV-type regularisation.

    This is an implementation of the classical PDHG method by Zhu
    and Chan [ZC2008]_. It aims to minimize

        ``1/2 * ||T(f) - g||^2 + lambda * ||S(f)||``

    where the first norm (data fit) is a Hilbert space norm and the
    second one (regularization) an arbitrary vector field norm.
    See Notes for mathematical details.

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
    dual : ``spars_op.range`` `element`
        Variable used in the dual iteration. It is updated in-place in
        each iteration step.
    spars_op : linear `Operator`
        Sparsifying operator ``S`` used in the regularization term.
        A typical choice is ``Gradient(primal.space)``.
    proj_op : `Operator`
        Projection operator applied to the dual in each step. A proper
        implementation of this step is crucial since it is the only
        place where the norm of the gradient enters the iteration.
        The projection operator must map the ``dual.space`` to
        itself.
    pstep, dstep : {`float`, `iterator`, `StepLength`}
        Rule for the selection of the primal/dual step parameter
        ``theta/tau``.
        A float is interpreted as constant step, an iterator is queried
        in each iteration for a new value, and a `StepLength` object
        can be used for an adaptive rule.
    reg_param : nonnegative `float`
        Regularisation parameter ``lambda`` steering the weight of the
        total variation penalizer. Larger value means less trust in the
        data and more weight on the regularizer.
    niter : positive `int`, optional
        Number of iterations to compute
    callback : `callable`, optional
        Function or object to be evaluated on the primal in each
        iteration. It is used as ``callback(xi)``, where ``xi`` is the
        current iterate.

    See also
    --------
    StepLength : rules for adaptive step length computation

    Notes
    -----

    The PDHG method solves for given data :math:`g\in\mathcal{Y}`,
    a linear forward operator
    :math:`\mathcal{T}:\mathcal{X}\\to\mathcal{Y}`, a sparsifying
    linear operator :math:`\mathcal{S}:\mathcal{X}\\to\mathcal{Z}`
    and a regularisation parameter :math:`\lambda > 0` the problem

        :math:`\min_{f \\in \mathcal{X}} \Big(\\frac{1}{2}
        \\lVert \mathcal{T}(f) - g\\rVert_{\mathcal{Y}}^2 +
        \\lambda \\lVert \mathcal{S}(f)\\rVert_{\mathcal{Z}} \Big)`,

    where the first norm is the Hilbert space norm in the data space
    and the second a norm in the range :math:`\mathcal{Z}` of the
    sparsifying opertor. For classical TV regularization, one can
    choose :math:`S = \\nabla` and
    :math:`\\lVert \cdot \\rVert_{\mathcal{Z}} = \\lVert \cdot \\rVert_1`.

    We assume that the norm in the regularising term can be written as

        :math:`\\lVert f\\rVert_{\mathcal{Z}} =
        \\sup_{\\lVert h\\rVert_\\ast \\leq 1}
        \\langle f, h\\rangle_{\mathcal{X}}`,

    with an *associate norm* :math:`\\lVert \cdot \\rVert_\\ast`.
    Then the method can be written as the following iteration scheme in
    the primal variable :math:`p` and the dual variable
    :math:`\mathbf{d}`:

        :math:`\mathbf{d}_{k+1} =
        \mathcal{P}_{\\lVert\cdot\\rVert_\\ast \\leq 1}
        \\big(\mathbf{d}_k - \\tau\, \mathcal{S}(f)\\big)`,

        :math:`p_{k+1} =
        p_k - \\theta \\big(\mathcal{T}^*\\big(\mathcal{T}(p_k)
        - g\\big) - \\lambda\, \mathcal{S}^*(\mathbf{d}_{k+1}) \\big)`.

    Here, :math:`\mathcal{P}_{\\lVert\cdot\\rVert_\\ast \\leq 1}` denotes
    the projection onto the unit ball in the associate norm.
    """
    # --- Handle input parameters --- #

    if not isinstance(fwd_op, Operator):
        raise TypeError('forward operator {!r} is not an Operator instance.'
                        ''.format(fwd_op))
    if not fwd_op.is_linear:
        raise ValueError('forward operator is not linear.')

    data = fwd_op.range.element(data)

    if primal not in fwd_op.domain:
        raise TypeError('primal {!r} is not an element of the forward '
                        'operator domain {!r}.'
                        ''.format(primal, fwd_op.domain))

    if not isinstance(spars_op, Operator):
        raise TypeError('sparsifying operator {!r} is not an Operator '
                        'instance.'.format(spars_op))
    if spars_op.domain != primal.space:
        raise OpDomainError('sparsifying operator domain is not the primal '
                            'variable space {!r}.'.format(primal.space))
    if dual not in spars_op.range:
        raise TypeError('dual {!r} is not an element of the range of the '
                        'sparsifying operator {}.'
                        ''.format(dual, spars_op.range))

    if not isinstance(proj_op, Operator):
        raise TypeError('projection operator {!r} is not an Operator instance.'
                        ''.format(proj_op))
    if proj_op.domain != dual.space:
        raise OpDomainError('projection operator domain is not the dual '
                            'variable space {!r}.'.format(dual.space))
    if proj_op.range != dual.space:
        raise OpRangeError('projection operator range is not the dual '
                           'variable space {!r}.'.format(dual.space))

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

    callback = kwargs.pop('callback', None)
    if callback is not None:
        if not callable(callback):
            raise TypeError('callback {!r} is not callable.'.format(callback))
        use_callback = True
    else:
        use_callback = False

    # --- Create auxiliary objects --- #

    res_op = ResidualOperator(fwd_op, data)
    tmp_primal = primal.space.element()
    tmp_dual = dual.space.element()
    tmp_data = fwd_op.range.element()

    # --- Run the method --- #

    for _ in range(niter):
        # Dual step
        # TODO: define dstep interface, currently ignores input
        cur_dstep = dstep(dual)
        spars_op(primal, out=tmp_dual)

        tmp_dual.lincomb(1, dual, -cur_dstep, tmp_dual)
        proj_op(tmp_dual, out=dual)

        # Primal step
        cur_pstep = pstep(primal)
        res_op(primal, out=tmp_data)
        # Doing the evaluation of the update parts for the primal in two
        # steps to save one temporary
        # TODO: extension to non-linear operators?
        fwd_op.adjoint(tmp_data, out=tmp_primal)
        primal.lincomb(1, primal, -cur_pstep, tmp_primal)
        spars_op.adjoint(dual, out=tmp_primal)
        primal.lincomb(1, primal, cur_pstep * reg_param, tmp_primal)

        if use_callback:
            callback(primal)

if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
