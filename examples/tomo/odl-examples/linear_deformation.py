# -*- coding: utf-8 -*-
"""
Example of creating an operator that acts as a linear transformation.
"""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

import odl
import numpy as np


def K(x, y, sigma):
    # Define the K matrix as symmetric gaussian
    return np.exp(-((x[0] - y[0])**2 + (x[1] - y[1])**2) / sigma**2) * np.eye(2)


def v(x, grid, alphas, sigma):
    """ Calculate the translation at point x """
    alpha1, alpha2 = alphas  # unpack translations per direction
    result = np.zeros_like(x)
    for i, (point, a1, a2) in enumerate(zip(grid.points(), alpha1, alpha2)):
        result += K(x, point, sigma).dot([a1, a2]).squeeze()
    return result


class LinearDeformation(odl.Operator):
    """ A linear deformation given by:
        ``the inverse of deformation \phi^-1(y) = y + v(y)``
        ``\phi.f(y) = f(\phi^-1(y))``

    Where ``f(x)`` is the input template, \phi.f(y) is the deformed template, 
    and ``v(y)`` is the inverse translation at point ``y``, 
    ``v(y)`` is computed using gaussian kernels with midpoints at ``grid``.
    
    The inverse translation ``v(y)`` is calculated as::
    
        v(y) = sum_k ( K(y, y_k) alpha_k )
    
    with ``alpha_k`` being a 2-dimensional vector for 2D case.
    
    Parameters
    ----------
    fspace : `DiscreteLp`
        Space of functions to be deformed
    vspace : `ProductSpace`
        Space of deformation coefficients ``alpha_k``. It has
        ``fspace.ndim`` components, each of the same dimension as
        ``fspace``.
    grid : `RegularGrid`
        Grid of the deformation control points ``x_k``. Its shape is equal
        to the shape of the components of ``vspace``.
    sigma : `float`
        Parameter for the Gaussian kernel
    """
    def __init__(self, fspace, vspace, grid, sigma, data):
        self.grid = grid
        self.sigma = sigma
        self.data = data
        super().__init__(domain=odl.ProductSpace(fspace, vspace),
                         range=fspace, linear=False)
    
    #Computing the deformed template: f(y + v(y))
    def _call(self, x):
        # Unpack input
        f, alphas = x
        # interpolant = f.space.extension(f.ntuple)  # this syntax is improved in pull #276

        # Array of output values
        out_values = np.zeros(f.size)

        for i, point in enumerate(self.range.points()):
            # Calculate deformation in each point
            point += v(point, self.grid, alphas, self.sigma)

            # Assign value at shifted point by interpolating. Set zero
            # if the shifted point is outside the rectangle.
            if point in self.domain[0].domain:
                # Evaluate interpolant at point
                out_values[i] = f.extension(point)
            else:
                # Zero boundary condition
                out_values[i] = 0

        return out_values
    
    #Computing the gradient of the shape regularization term: 2*K*alpha
    def grad_shape_regular_term(self, x):
        f, alphas = x
        out_values = np.zeros(f.size)
        for i, point in enumerate(self.range.points()):
            out_values= 2*v(point, self.grid, alphas, self.sigma)
        return out_values
        
    #Computing the gradient of the regularization term on deformed template: \nabla R(f())      
    def grad_image_regular_term(self, x):
        pass
    
    #Computing the gradient of the fitting term
    def grad_fittting_term(self, x):
        pass
        
    

# Discretization of the space
m = 101  # Number of gridpoints for discretization
spc = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [m, m], dtype='float32')

# deformation space
n = 5  # number of gridpoints for deformation, usually smaller than m
vspace = odl.ProductSpace(odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [n, n]), 2)

# Deformation operator
deformation = LinearDeformation(spc, vspace, vspace[0].grid, sigma=0.2, data=0)

# Create input function
f = odl.util.shepp_logan(spc, True)

# Create data
detector_partition = odl.uniform_partition(-0.75, 0.75, 151)

angle_interval = odl.Interval(0, np.pi)
angle_grid = odl.TensorGrid([0, np.pi/4, np.pi/2])
angle_partition = odl.RectPartition(angle_interval, angle_grid)

#angle_partition = odl.uniform_partition(0, 2 * np.pi, 3)

geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)

xray_trafo = odl.tomo.RayTransform(spc, geometry, impl='astra_cuda')

proj_data = xray_trafo(f)
backproj = xray_trafo.adjoint(proj_data)

# Shows a slice of the phantom, projections, and reconstruction
proj_data.show(title='Projection data (sinogram)')
backproj.show(title='Back-projected data')


# Create deformation field
values = np.zeros([2, n, n])
values[0, :, :n//2] = 0.1  # movement in "x" direction
values[1, n//2, :] = 0.05   # movement in "y" direction
def_coeff = vspace.element(values)

# Show input
f.show(title='f')
def_coeff.show(title='deformation')

# Calculate deformed function
result = deformation.range.element()
deformation([f, def_coeff], out=result)
result.show(title='result')



