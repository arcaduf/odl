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

"""Denoising using the Chambolle-Pock solver with TV & entropy-type data term.

Solves the optimization problem with the Kullback-Leibler data divergence

    min_{x > 0}  sum(A(x) - g ln(A(x)) + lam || |grad(x)| ||_1

For details see :ref:`chambolle_pock`, :ref:`proximal_operators`, and
references therein.
"""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

import numpy as np
import scipy
import scipy.ndimage

import odl


# Read test image:
# convert integer values to float, and rotate to get the image upright
image = np.rot90(scipy.misc.ascent()[::2, ::2], 3).astype('float')
shape = image.shape

# Rescale
image *= 100 / image.max()

# Add noise
noisy_image = np.random.poisson(1 + image)

# Discretized spaces and vectors
discr_space = odl.uniform_discr([0, 0], shape, shape)
orig = discr_space.element(image)
noisy = discr_space.element(noisy_image)


# --- Set up the inverse problem --- #


# Gradient operator
gradient = odl.Gradient(discr_space, method='forward')

# Matrix of operators
op = odl.BroadcastOperator(odl.IdentityOperator(discr_space), gradient)


# Proximal operator related to the primal variable

# Non-negativity constraint
proximal_primal = odl.solvers.proximal_nonnegativity(op.domain)

# Proximal operators related to the dual variable

# l2-data matching
prox_convconj_kl = odl.solvers.proximal_convexconjugate_kl(discr_space,
                                                           lam=1.0, g=noisy)

# TV-regularization: l1-semi norm of grad(x)
prox_convconj_l1 = odl.solvers.proximal_convexconjugate_l1(gradient.range,
                                                           lam=0.1)

# Combine proximal operators: the order must match the order of operators in K
proximal_dual = odl.solvers.combine_proximals([prox_convconj_kl,
                                               prox_convconj_l1])

# Optional: pass partial objects to solver
partial = (odl.solvers.PrintIterationPartial() &
           odl.solvers.ShowPartial())


# --- Select solver parameters and solve using Chambolle-Pock --- #


# Estimated operator norm, add 10 percent to ensure ||K||_2^2 * sigma * tau < 1
op_norm = 1.1 * odl.operator.oputils.power_method_opnorm(op, 100)

niter = 100  # Number of iterations
tau = 10.0 / op_norm  # Step size for the primal variable
sigma = 0.1 / op_norm  # Step size for the dual variable

# Starting point
x = op.domain.zero()

# Run algorithms (and display intermediates)
odl.solvers.chambolle_pock_solver(
    op, x, tau=tau, sigma=sigma, proximal_primal=proximal_primal,
    proximal_dual=proximal_dual, niter=niter, partial=partial)

# Display images
orig.show(title='original image')
noisy.show(title='noisy image')
x.show(title='denoised', show=True)  # show and hold
