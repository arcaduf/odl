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

"""Denoising using the PDHG solver with TV regularization term.

Solves the optimization problem with the Kullback-Leibler data divergence

    min_{x > 0}  sum(A(x) - g ln(A(x)) + lam || |grad(x)| ||_1

For details see :ref:`chambolle_pock`, :ref:`proximal_operators`, and
references therein.
"""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

import numpy as np
import scipy
import scipy.ndimage

import odl

# convert integer values to float, and rotate to get the image upright
image = np.rot90(scipy.misc.ascent()[::2, ::2], 3).astype('float')
image = np.rot90(scipy.misc.face(gray=True)[::-1, ::-1]).astype('float')
#reco_space = odl.uniform_discr([-1, -1], [1, 1], (512, 512))
#image = odl.util.phantom.shepp_logan(reco_space, modified=True)

# Rescale
image *= 100 / np.amax(image)

# Add noise
noisy_image = image + np.random.normal(0, scale=10.0, size=image.shape)

# Discretized spaces and vectors
reco_space = odl.uniform_discr([0, 0], image.shape, image.shape)
orig = reco_space.element(image)
noisy = reco_space.element(noisy_image)


# --- Set up the inverse problem --- #


# Forward operator is the identity
fwd_op = odl.IdentityOperator(reco_space)

# Sparsifying operator - the spatial gradient
spars_op = odl.Gradient(reco_space, method='forward')


# Projection operator
class NormalizationOperator(odl.Operator):
    """Operator normalizing a vector field with its 1-norm."""
    def __init__(self, space, exponent):
        assert isinstance(space, odl.ProductSpace)
        assert isinstance(space[0], odl.DiscreteLp)
        super().__init__(domain=space, range=space, linear=False)
        self.exponent = float(exponent)

    def _call(self, x):
        exponent = odl.util.utility.conj_exponent(self.exponent)
        const = x[0].cell_volume if np.isfinite(exponent) else 1.0
        weight = odl.space.ntuples.FnConstWeighting(const,
                                                    exponent=exponent)
        pw_norm = odl.PointwiseNorm(self.domain)
        x_norm = weight.norm(pw_norm(x).ntuple)
        if x_norm <= 1.0:
            return x
        else:
            return x / x_norm


class NormProjectionOperator(odl.discr.tensor_ops.PointwiseTensorFieldOperator):

    """Projection onto ``||.|| <= 1``."""

    def __init__(self, space, exponent, proj_iter=2, descent_method=None):
        if isinstance(space, odl.ProductSpace):
            self.pw_norm = odl.PointwiseNorm(space)
        elif isinstance(space, odl.DiscreteLp):
            self.pw_norm = None
        else:
            raise TypeError

        super().__init__(domain=space, range=space, linear=True)

        # TODO: adapt for variable exponent
        self.exponent = float(exponent)
        self.conj_exponent = odl.util.utility.conj_exponent(self.exponent)
        if np.isfinite(self.conj_exponent):
            const = self.base_space.cell_volume
        else:
            const = 1.0
        self.fn_norm = odl.weighted_norm(const, exponent=self.conj_exponent)
        self.proj_iter = int(proj_iter)
        if descent_method is None:
            self.descent_method = odl.solvers.bfgs_method
        else:
            if not callable(descent_method):
                raise TypeError
            self.descent_method = descent_method

    def _conj_exp_norm(self, x):
        if self.pw_norm is not None:
            y = self.pw_norm(x)
        else:
            y = x
        # TODO: boundary behavior is not respected. Better to define a new
        # space with the conjugate exponent and use its norm.
        return self.fn_norm(y.ntuple)

    def _proj_inf_norm(self, x):
        # TODO: Shrinkage
        raise NotImplementedError

    def _proj_2_norm(self, x):
        # TODO: Should be simple?
        raise NotImplementedError

    def _call(self, x, out):
        def cost(y):
            y_norm = self._conj_exp_norm(y)
            if y_norm > 1:
                return np.inf
            else:
                return x.dist(y)

        line_search = odl.solvers.BacktrackingLineSearch(cost)
        cost_grad = odl.ResidualOperator(
            odl.IdentityOperator(self.domain), x)

        out[:] = self.domain.zero()
        self.descent_method(cost_grad, out, line_search, niter=self.proj_iter)


p = 1.0
#proj_op = NormalizationOperator(spars_op.range, exponent=p)
proj_op = NormProjectionOperator(spars_op.range, exponent=p)

# Initialize the primal and dual variables. We need to choose something
# that does not yield Gradient(primal) = 0 and dual = 0.
primal = odl.util.phantom.submarine_phantom(spars_op.domain)
dual = spars_op.range.zero()


# --- Run the algorithm --- #

partial = odl.solvers.ShowPartial()

#odl.solvers.primal_dual_hybrid_gradient(
#    fwd_op=odl.IdentityOperator(discr_space), data=noisy, primal=primal,
#    dual=dual, spars_op=spars_op, proj_op=proj_op,
#    pstep=0.5, dstep=0.000001, reg_param=400000, niter=20, callback=partial)

reg_param = spars_op.domain.cell_volume * 0.5
odl.solvers.primal_dual_hybrid_gradient(
    fwd_op=fwd_op, data=noisy, primal=primal, dual=dual, spars_op=spars_op,
    proj_op=proj_op, pstep=0.2, dstep=0.1, reg_param=reg_param, niter=20,
    callback=partial, balance=True)

# Display images
orig.show(title='original image')
noisy.show(title='noisy image')
primal.show(title='denoised, p={}'.format(p), show=True)  # show and hold
