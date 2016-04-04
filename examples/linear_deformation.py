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

"""Example collection for linearized deformations."""

from __future__ import print_function, division
from functools import partial
import numpy as np

import odl


# kernel function for any dimensional
def gauss_kernel(x, sigma):
    return np.exp(-x ** 2 / (2 * sigma ** 2))

# Fix the sigma parameter in the kernel
sigma = 0.3

# Discretization of the space
m = 101  # Number of gridpoints for discretization
spc = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [m, m], dtype='float32')

# deformation space
n = 5  # number of gridpoints for deformation, usually smaller than m
cptssapce = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [n, n])
vspace = odl.ProductSpace(cptssapce, 2)

# Deformation operator
#deformation = LinearDeformation(spc, vspace, vspace[0].grid, sigma=0.2, data=0)

# Create input function as Shepp-Logan phantom
template = odl.util.shepp_logan(spc, True)

# Create target function as submarine phantom
target = odl.util.submarine_phantom(spc, smooth=True, taper=50.0)

# Create projection domain
detector_partition = odl.uniform_partition(-0.75, 0.75, 151)

# Create projection directions
angle_interval = odl.Interval(0, np.pi)
angle_grid = odl.TensorGrid([0, np.pi/4, np.pi/2])
angle_partition = odl.RectPartition(angle_interval, angle_grid)

# Create uniform direction in [0, 2*pi]
# angle_partition = odl.uniform_partition(0, 2 * np.pi, 3)

# Create 2D parallel projection geometry
geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)

# Create forward projections by Radon transform
xray_trafo_op = odl.tomo.RayTransform(spc, geometry, impl='astra_cuda')
proj_data = xray_trafo_op(template)

# Create backprojection
backproj = xray_trafo_op.adjoint(proj_data)

# Shows a slice of the phantom, projections, and reconstruction
proj_data.show(title='Projection data (sinogram)')
backproj.show(title='Back-projected data')

# Create deformation field
values = np.zeros([2, n, n])
values[0, :, :n//2] = 0.025  # movement in "x" direction
values[1, n//2, :] = 0.01   # movement in "y" direction
def_coeff = vspace.element(values)

# Show input
template.show(title='Template')
def_coeff.show(title='deformation')

kernel = partial(gauss_kernel, sigma=sigma)
tpl_op = TemplateDeformationOperator(vspace, cptssapce.grid, template, kernel)
deformed_image = tpl_op(def_coeff) # deformed image
deformed_image.show(title='deformed image')

tpl_op_deriv = tpl_op.derivative(def_coeff)  # an operator
beta = np.ones_like(def_coeff)
deriv_result = tpl_op_deriv(beta)  # the same betas as alphas
deriv_result.show('derivative result')


#######################test 1##########################
# Compute the gradient of L2 fitting term by mathematical method directly
proj_data_deformed_image = xray_trafo_op(deformed_image)
backproj_diff = xray_trafo_op.adjoint(proj_data_deformed_image - proj_data)
backproj_diff.show('back projection of difference')

tpl_deriv_adj = tpl_op_deriv.adjoint
adj_result = tpl_deriv_adj(backproj_diff) # Back projection of difference as the function f(x)
adj_result.show('adjoint result')


########################test 2#########################
# Composition of the L2 fitting term with a deformation operator
l2_data_fit = L2DataMatchingFunctional(xray_trafo_op.range, proj_data)
data_fitting_term = l2_data_fit * xray_trafo_op * tpl_op

grad = data_fitting_term.gradient(def_coeff)
grad.show('gradient of fitting term')

## compare adj_result with grad
#(adj_result - grad).show('difference')

# compute the gradient of shape-based regularization term
kernelmatrix = gaussian_kernel_matrix(cptssapce.grid, sigma)
shape_func = ShapeRegularizationFunctional(vspace, kernelmatrix)
grad_shape_func = shape_func._gradient(def_coeff)
grad_shape_func.show('gradient of shape functional')






##################################################################################################################
#################################################simple test######################################################
## quite simple testing
#
#spc = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [2, 2], dtype='float32')
#template = spc.element([0, 1, 1, 0])
#template.show('Template')
#cptssapce = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [1, 1])
#vspace = odl.ProductSpace(cptssapce, 2)
#def_coeff = vspace.element([0, 1])
#sigma = 0.3
#kernel = partial(gauss_kernel, sigma=sigma)
#tpl_op = TemplateDeformationOperator(vspace, cptssapce.grid, template, kernel)
#deformed_image = tpl_op(def_coeff) # deformed image
#deformed_image.show('deformed image')
#
#kernelmatrix = gaussian_kernel_matrix(cptssapce.grid, sigma)
#shape_func = ShapeRegularizationFunctional(vspace, kernelmatrix)
#grad_shape_func = shape_func._gradient(def_coeff)
##grad_shape_func.show('gradient of shape functional')
#
#
### test code
###(tpl_op(vspace.one() + 0.1 * vspace.one()) - tpl_op(vspace.one()) - 0.1 * tpl_op.derivative(vspace.one())(vspace.one())).norm()

###########################################################################################################


## Calculate deformed function
#result = deformation.range.element()
#deformation([template, def_coeff], out=result)
#result.show(title='result')


#def _deform(self, template, alphas):
#        # Array of output values
#        out_values = np.zeros(template.size)
#
#        for i, point in enumerate(self.template.space.points()):
#            # Calculate deformation in each point
#            point += v(point, self.grid, alphas, self.sigma)
#
#            # Assign value at shifted point by interpolating. Set zero
#            # if the shifted point is outside the rectangle.
#            if point in self.domain[0].domain:
#                # Evaluate interpolant at point
#                out_values[i] = template.extension(point)
#            else:
#                # Zero boundary condition
#                out_values[i] = 0
#
#        return out_values
