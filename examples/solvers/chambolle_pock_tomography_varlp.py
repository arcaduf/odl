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

"""Total variation tomography using the Chambolle-Pock solver.

Solves the optimization problem

    min_x  1/2 ||A(x) - g||_2^2 + lam || |grad(x)| ||_1

For details see :ref:`chambolle_pock`, :ref:`proximal_operators`, and
references therein.
"""

import numpy as np
import odl


reco_space = odl.uniform_discr([-10, -10], [10, 10], (300, 300),
                               dtype='float32')
#small_space = odl.uniform_discr([-10, -10], [10, 10], (100, 100),
#                                dtype='float32')
#tmp = odl.util.phantom.submarine_phantom(small_space, smooth=False)
#tmp += odl.util.phantom.submarine_phantom(small_space, smooth=True, taper=5)
#phantom = reco_space.zero()
#phantom.asarray()[100:200, 100:200] = tmp
phantom = odl.util.phantom.submarine_phantom(reco_space, smooth=False)
phantom += odl.util.phantom.submarine_phantom(reco_space, smooth=True, taper=5)
phantom.show('Phantom')

# Define the exponent: we use a convolution with a 3x3 constant kernel to
# broaden the boundaries (safety margin) without tapering off. After
# that, we set areas below a threshold to 2.0, above to 1.0.
# Possible issues:
# - Small jumps have smaller Laplacian -> low threshold needed, otherwise miss
# - Threshold is a parameter. How to coose?
# - Sensitivity to noise
# - Binary image - perhaps better to have a smooth function

s = 0.5


def exp_kernel(x):
    scaled = [xi / (np.sqrt(2) * s) for xi in x]
    return np.exp(-sum(xi ** 2 for xi in scaled))


def p(x):
    return (((np.abs(x[0]) < 5.2) & (np.abs(x[1]) < 5.2)).astype(float) -
            ((np.abs(x[0]) < 4.8) & (np.abs(x[1]) < 4.8)))


add_kernel = np.ones((10, 10))
exp_conv = odl.Convolution(reco_space, exp_kernel, impl='default_ft',
                           scale=False)
add_conv = odl.Convolution(reco_space, add_kernel, impl='scipy_convolve',
                           scale=False)
lapl = odl.Laplacian(reco_space)
abs_lapl = np.abs(lapl(phantom))
conv_abs_lapl = np.maximum(exp_conv(abs_lapl), 0)
#conv_abs_lapl[:] = add_conv(conv_abs_lapl)
conv_abs_lapl *= 1.5 / np.max(conv_abs_lapl)
conv_abs_lapl[:] = np.minimum(conv_abs_lapl, 1)
var_exponent = 2.0 - conv_abs_lapl
var_exponent_arr = var_exponent.asarray()

q = reco_space.element(p)
conv_q = exp_conv(q)
conv_q *= 1.5 / np.max(conv_q)
conv_q[:] = np.minimum(conv_q, 1)
var_exponent = 2.0 - conv_q


# Shift the exponent function (disturbance)
new_exp = 2 * np.ones_like(var_exponent)
new_exp[:, 10:] = var_exponent_arr[:, :-10]
var_exponent[:] = new_exp

# Disturb the exponent function
#var_exponent_arr[:, 100:150] = 2.0

var_exponent.show('Exponent function')


# --- Set up the forward operator (ray transform) --- #


# Make a parallel beam geometry with flat detector
# Angles: uniformly spaced, n = 360, min = 0, max = 2 * pi
angle_partition = odl.uniform_partition(0, 2 * np.pi, 360)
# Detector: uniformly sampled, n = 558, min = -30, max = 30
detector_partition = odl.uniform_partition(-30, 30, 558)
geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)

# The implementation of the ray transform to use, options:
# 'scikit'                    Requires scikit-image (can be installed by
#                             running ``pip install scikit-image``).
# 'astra_cpu', 'astra_cuda'   Require astra tomography to be installed.
#                             Astra is much faster than scikit. Webpage:
#                             https://github.com/astra-toolbox/astra-toolbox
impl = 'astra_cuda'

# Ray transform as forward operator
ray_trafo = odl.tomo.RayTransform(reco_space, geometry, impl=impl)


# --- Generate data --- #

data = ray_trafo(phantom)
noisy_data = data + 3 * odl.util.phantom.white_noise(data.space)


# --- Set up the inverse problem --- #


# Initialize gradient operator
gradient = odl.Gradient(reco_space, method='forward')

# Column vector of two operators
op = odl.BroadcastOperator(ray_trafo, gradient)

# Create the proximal operator for unconstrained primal variable
proximal_primal = odl.solvers.proximal_zero(op.domain)

# Create proximal operators for the dual variable

# L2-data matching
#prox_convconj_l2_data = odl.solvers.proximal_convexconjugate_l2(
#    ray_trafo.range, g=data)

prox_convconj_l2_data = odl.solvers.proximal_convexconjugate_l2_squared(
    ray_trafo.range, g=noisy_data)

# TV-regularization with variable Lp
lam = 0.1
prox_var_lp = odl.solvers.proximal_variable_lp(gradient.range, var_exponent,
                                               lam=lam)
prox_convconj_varlp = odl.solvers.proximal_convexconjugate(prox_var_lp)
prox_convconj_l1 = odl.solvers.proximal_convexconjugate_l1(
    gradient.range, lam=lam)

prox_convconj_l2 = odl.solvers.proximal_convexconjugate_l2_squared(
    gradient.range, lam=lam)

# Combine proximal operators, order must correspond to the operator K
proximal_dual = odl.solvers.combine_proximals(
    [prox_convconj_l2_data, prox_convconj_varlp])
#proximal_dual = odl.solvers.combine_proximals(
#    [prox_convconj_l2_data, prox_convconj_l2])


# --- Select solver parameters and solve using Chambolle-Pock --- #


# Estimated operator norm, add 10 percent to ensure ||K||_2^2 * sigma * tau < 1
op_norm = 1.5 * odl.operator.oputils.power_method_opnorm(op, 5)
print('operator norm estimate: ', op_norm)

niter = 2000  # Number of iterations
tau = 1.0  # Step size for the primal variable
sigma = 1.0 / (op_norm ** 2 * tau)  # Step size for the dual variable
gamma = 0.1
theta = 0.85

# Optionally pass partial to the solver to display intermediate results
partial = (odl.solvers.util.PrintIterationPartial() &
           odl.solvers.util.ShowPartial(display_step=50))

# Choose a starting point
x = op.domain.zero()

# Run the algorithm
odl.solvers.chambolle_pock_solver(
    op, x, tau=tau, sigma=sigma, gamma=gamma, theta=theta,
    proximal_primal=proximal_primal,
    proximal_dual=proximal_dual, niter=niter, partial=partial)

# Display images
phantom.show(title='Phantom')
noisy_data.show(title='Data')
x.show(title='Reconstruction', show=True)
