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

"""Total variation denoising using the Chambolle-Pock solver.

Solves the optimization problem

    min_{x >= 0}  1/2 ||x - g||_2^2 + lam || |grad(x)| ||_1

Where ``grad`` the spatial gradient and ``g`` is given noisy data.

For further details and a description of the solution method used, see
:ref:`chambolle_pock` in the ODL documentation.
"""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

import numpy as np
import matplotlib.pyplot as plt
import odl

# --- Set up the forward operator (ray transform) --- #


# Discrete reconstruction space: discretized functions on the rectangle
# [-20, 20]^2 with 300 samples per dimension.
reco_space = odl.uniform_discr(
    min_corner=[-20, -20], max_corner=[20, 20], nsamples=[300, 300],
    dtype='float32')

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

# Ray transform aka forward projection.
ray_trafo = odl.tomo.RayTransform(reco_space, geometry, impl=impl)


# --- Generate artificial data --- #


# Create phantom
phantom = odl.util.shepp_logan(reco_space, modified=True)

# Create sinogram of forward projected phantom with noise
data = ray_trafo(phantom)
data += odl.util.white_noise(ray_trafo.range) * np.mean(data) * 0.02


# --- Set up the inverse problem --- #


# Initialize gradient operator
gradient = odl.Gradient(reco_space, method='backward')

# Column vector of two operators
op = odl.BroadcastOperator(ray_trafo, gradient)

# Create the proximal operator for unconstrained primal variable
proximal_primal = odl.solvers.proximal_zero(op.domain)

# Create proximal operators for the dual variable

# l2-data matching
prox_convconj_l2 = odl.solvers.proximal_convexconjugate_l2(ray_trafo.range,
                                                           g=data)

# TV-regularization i.e. the l1-norm
prox_convconj_l1 = odl.solvers.proximal_convexconjugate_l1(
    gradient.range, lam=0.01)

# Combine proximal operators, order must correspond to the operator K
proximal_dual = odl.solvers.combine_proximals(
    [prox_convconj_l2, prox_convconj_l1])


# Set some general parameters
op_norm_ray_trafo = 1.0 * odl.power_method_opnorm(ray_trafo, 10, phantom)
op_norm_gradient = 1.0 * odl.power_method_opnorm(gradient, 100, phantom)
print('Operator norms, Ray trafo: {}, gradient: {}'.format(op_norm_ray_trafo,
                                                           op_norm_gradient))

# --- Run algorithms without preconditioner

niter = 100

# Create a function to save the partial errors
partial = odl.solvers.StorePartial(function=lambda x: (x-phantom).norm())
#partial &= odl.solvers.ShowPartial()

# Step sizes
op_norm = max(op_norm_ray_trafo, op_norm_gradient)
tau = sigma = 1.0 / op_norm

x = op.domain.zero()  # Starting point
odl.solvers.chambolle_pock_solver(
    op, x, tau=tau, sigma=sigma, proximal_primal=proximal_primal,
    proximal_dual=proximal_dual, niter=niter, partial=partial)


# --- Run algorithm with preconditoner

# Fourier transform in detector direction
fourier = odl.trafos.FourierTransform(ray_trafo.range, axes=[1])

# Create ramp in the detector direction
cut = 25.0
c = np.pi / cut
ramp_function = fourier.range.element(lambda x: 0*x[0] + 5 + np.abs(np.sin(x[1] * c) / c) * (np.abs(x[1]) < cut))

# Create ramp filter via the convolution formula with fourier transforms
ramp_filter = (fourier.inverse * ramp_function * fourier)

filtered_ray_trafo = ramp_filter * ray_trafo

# Column vector of two operators
op = odl.BroadcastOperator(filtered_ray_trafo, gradient)

# Create proximal operators for the dual variable

# l2-data matching
prox_convconj_l2 = odl.solvers.proximal_convexconjugate_l2(filtered_ray_trafo.range,
                                                           g=ramp_filter(data))

# TV-regularization i.e. the l1-norm
prox_convconj_l1 = odl.solvers.proximal_convexconjugate_l1(
    gradient.range, lam=0.01)

# Combine proximal operators, order must correspond to the operator K
proximal_dual = odl.solvers.combine_proximals(
    [prox_convconj_l2, prox_convconj_l1])

# Create a function to save the partial errors
partial_prec = odl.solvers.StorePartial(function=lambda x: (x-phantom).norm())
partial_prec_show = partial_prec & odl.solvers.ShowPartial()

# Set some general parameters

# Step sizes
op_norm_ray_trafo = 1.0 * odl.power_method_opnorm(filtered_ray_trafo, 10, phantom)
op_norm_gradient = 1.0 * odl.power_method_opnorm(gradient, 100, phantom)
print('Operator norms, F Ray trafo: {}, gradient: {}'.format(op_norm_ray_trafo,
                                                             op_norm_gradient))
op_norm = max(op_norm_ray_trafo, op_norm_gradient)
tau = sigma = 1.0 / op_norm


x_precon = op.domain.zero()  # Starting point
odl.solvers.chambolle_pock_solver(
    op, x_precon, tau=tau, sigma=sigma, proximal_primal=proximal_primal,
    proximal_dual=proximal_dual, niter=niter, partial=partial_prec_show)

# results
x.show('Standard')
x_precon.show('With preconditioner')

plt.figure()
plt.loglog(np.arange(niter), partial, label='Standard')
plt.loglog(np.arange(niter), partial_prec, label='With preconditioner')
plt.legend(loc=3)
plt.xlabel('Iteration')
plt.ylabel('Error norm')
plt.title('Convergence')
