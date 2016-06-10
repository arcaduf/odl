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

"""
Example of shape-based image reconstruction using linearized deformations.
"""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
from builtins import super
from numbers import Number
from odl.operator.operator import Operator, OperatorComp
import odl
import numpy as np
import time
import matplotlib.pyplot as plt
standard_library.install_aliases()


class Functional(Operator):

    """Quick hack for a functional class."""

    def __init__(self, domain, linear=False):
        """Initialize a new instance.

        Parameters
        ----------
        domain : `LinearSpace`
            Set of elements on which the functional can be evaluated
        linear : `bool`, optional
            If `True`, assume that the functional is linear
        """
        super().__init__(domain=domain, range=domain.field, linear=linear)

    def gradient(self, x, out=None):
        """Evaluate the gradient of the functional.

        Parameters
        ----------
        x : domain element-like
            Point in which to evaluate the gradient
        out : domain element, optional
            Element into which the result is written

        Returns
        -------
        out : domain element
            Result of the gradient calcuation. If ``out`` was given,
            the returned object is a reference to it.
        """
        raise NotImplementedError

    def __mul__(self, other):
        """Return ``self * other``.

        If ``other`` is an operator, this corresponds to
        operator composition:

            ``op1 * op2 <==> (x --> op1(op2(x))``

        If ``other`` is a scalar, this corresponds to right
        multiplication of scalars with operators:

            ``op * scalar <==> (x --> op(scalar * x))``

        If ``other`` is a vector, this corresponds to right
        multiplication of vectors with operators:

            ``op * vector <==> (x --> op(vector * x))``

        Note that left and right multiplications are generally
        different.

        Parameters
        ----------
        other : {`Operator`, `LinearSpaceVector`, scalar}
            `Operator`:
                The `Operator.domain` of ``other`` must match this
                operator's `Operator.range`.

            `LinearSpaceVector`:
                ``other`` must be an element of this operator's
                `Operator.domain`.

            scalar:
                The `Operator.domain` of this operator must be a
                `LinearSpace` and ``other`` must be an
                element of the ``field`` of this operator's
                `Operator.domain`.

        Returns
        -------
        mul : `Functional`
            Multiplication result

            If ``other`` is an `Operator`, ``mul`` is a
            `FunctionalComp`.

            If ``other`` is a scalar, ``mul`` is a
            `FunctionalRightScalarMult`.

            If ``other`` is a vector, ``mul`` is a
            `FunctionalRightVectorMult`.

        """
        if isinstance(other, Operator):
            return FunctionalComp(self, other)
        elif isinstance(other, Number):
            # Left multiplication is more efficient, so we can use this in the
            # case of linear operator.
            raise NotImplementedError
            if self.is_linear:
                return OperatorLeftScalarMult(self, other)
            else:
                return OperatorRightScalarMult(self, other)
        elif isinstance(other, LinearSpaceVector) and other in self.domain:
            raise NotImplementedError
            return OperatorRightVectorMult(self, other.copy())
        else:
            return NotImplemented


class FunctionalComp(Functional, OperatorComp):

    """Composition of a functional with an operator."""

    def __init__(self, func, op, tmp1=None, tmp2=None):
        """Initialize a new instance.

        Parameters
        ----------
        func : `Functional`
            The left ("outer") operator
        op : `Operator`
            The right ("inner") operator. Its range must coincide with the
            domain of ``func``.
        tmp1 : `element` of the range of ``op``, optional
            Used to avoid the creation of a temporary when applying ``op``
        tmp2 : `element` of the range of ``op``, optional
            Used to avoid the creation of a temporary when applying the
            gradient of ``func``
        """
        if not isinstance(func, Functional):
            raise TypeError('functional {!r} is not a Functional instance.'
                            ''.format(func))

        OperatorComp.__init__(self, left=func, right=op, tmp=tmp1)

        if tmp2 is not None and tmp2 not in self._left.domain:
            raise TypeError('second temporary {!r} not in the domain '
                            '{!r} of the functional.'
                            ''.format(tmp2, self._left.domain))
        self._tmp2 = tmp2

    def gradient(self, x, out=None):
        """Gradient of the compositon according to the chain rule.

        Parameters
        ----------
        x : domain element-like
            Point in which to evaluate the gradient
        out : domain element, optional
            Element into which the result is written

        Returns
        -------
        out : domain element
            Result of the gradient calcuation. If ``out`` was given,
            the returned object is a reference to it.
        """
        if out is None:
            # adj_op = self._right.derivative(x).adjoint
            # return adj_op(self._left.gradient(self._right(x)))
            return self._right.derivative(x).adjoint(
                self._left.gradient(self._right(x)))
        else:
            if self._tmp is not None:
                tmp_op_ran = self._right(x, out=self._tmp)
            else:
                tmp_op_ran = self._right(x)

            if self._tmp2 is not None:
                tmp_dom = self._left.gradient(tmp_op_ran, out=self._tmp2)
            else:
                tmp_dom = self._left.gradient(tmp_op_ran)

            self._right.derivative(x).adjoint(tmp_dom, out=out)


class DisplacementOperator(Operator):

    """Operator mapping parameters to an inverse displacement for
    the domain of the target.

    This operator computes for the momenta::

        alpha --> D(alpha)

    where

        D(alpha)(y) = v(y).

    The operator vector field ``v`` depends on the deformation parameters
    ``alpha_j`` as follows::

        v(y) = sum_j (K(y, y_j) * alpha_j)

    Here, ``K`` is the RKHS kernel matrix and each ``alpha_j`` is an
    element of ``R^n``, can be seen as the momenta alpha at control point y_j.
    """

    def __init__(self, par_space, control_points, discr_space, ft_kernel):
        """Initialize a new instance.

        Parameters
        ----------
        par_space : `ProductSpace` or `Rn`
            Space of the parameters. For one-dimensional deformations,
            `Rn` can be used. Otherwise, a `ProductSpace` with ``n``
            components is expected.
        control_points : `TensorGrid` or `array-like`
            The points ``x_j`` controlling the deformation. They can
            be given either as a tensor grid or as a point array. In
            the latter case, its shape must be ``(N, n)``, where
            ``n`` is the dimension of the template space, and ``N``
            the number of ``alpha_j``, i.e. the size of (each
            component of) ``par_space``.
        discr_space : `DiscreteSpace`
            Space of the image grid of the template.
        kernel : `callable`
            Function to determine the kernel at the control points ``K(y_j)``
            The function must accept a real variable and return a real number.
        """
        # TODO: use kernel operator instead of function & matrix
        if isinstance(par_space, odl.Fn):
            # Make a product space with one component
            par_space = odl.ProductSpace(par_space, 1)
        elif isinstance(par_space, odl.ProductSpace):
            pass
        else:
            raise TypeError('expected Rn or ProductSpace as par_space, got '
                            '{!r}.'.format(par_space))

        if par_space.size != discr_space.ndim:
            raise ValueError('dimensions of product space and image grid space'
                             ' do not match ({} != {})'
                             ''.format(par_space.size, discr_space.ndim))

        # The operator maps from the parameter space to an inverse
        # displacement for the domain of the target.

        self.discr_space = discr_space
        self.range_space = odl.ProductSpace(self.discr_space,
                                            self.discr_space.ndim)

        super().__init__(par_space, self.range_space, linear=True)

        self.ft_kernel = ft_kernel

        if not isinstance(control_points, odl.TensorGrid):
            self._control_pts = np.asarray(control_points)
            if self._control_pts.shape != (self.num_contr_pts, self.ndim):
                raise ValueError(
                    'expected control point array of shape {}, got {}.'
                    ''.format((self.num_contr_pts, self.ndim),
                              self.control_points.shape))
        else:
            self._control_pts = control_points

        # TODO: check that number of control points is the same as alphas

    @property
    def ndim(self):
        """Number of dimensions of the deformation."""
        return self.domain.size

    @property
    def contr_pts_is_grid(self):
        """`True` if the control points are given as a grid."""
        return isinstance(self.control_points, odl.TensorGrid)

    @property
    def control_points(self):
        """The control points for the deformations."""
        return self._control_pts

    @property
    def num_contr_pts(self):
        """Number of control points for the deformations."""
        if self.contr_pts_is_grid:
            return self.control_points.size
        else:
            return len(self.control_points)

    @property
    def image_grid(self):
        """Spatial sampling grid of the image space."""
        return discr_space.grid

    def displacement_ft(self, alphas):
        """Calculate the inverse translation at point y by n-D FFT.

        inverse displacement: v(y)

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points. Note that here N = M.
        kernel: Gaussian kernel fuction
        """
        ft_momenta = vectorial_ft_fit_op(alphas)
        ft_displacement = self.ft_kernel * ft_momenta
        return vectorial_ft_fit_op.inverse(ft_displacement)
        # scaling
#        return (vectorial_ft_op_inverse(ft_displacement) /
#                self.discr_space.cell_volume * 2.0 * np.pi)

    def _call(self, alphas):
        """Implementation of ``self(alphas, out)``.

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.

        out : `DiscreteLpVector`
            Element where the result is stored
        """

        return self.displacement_ft(alphas)

    def derivative(self, alphas):
        """Frechet derivative of this operator in ``alphas``.

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.

        Returns
        -------
        deriv : `Operator`
            The derivative of this operator, evaluated at ``alphas``
        """
        deriv_op = DisplacementDerivative(
            alphas, self.control_points, self.discr_space, self.ft_kernel)

        return deriv_op


class DisplacementDerivative(DisplacementOperator):

    """Frechet derivative of the displacement operator at alphas."""

    def __init__(self, alphas, control_points, discr_space, ft_kernel):
        """Initialize a new instance.

        Parameters
        ----------
        alphas : `ProductSpaceVector`
            Displacement parameters in which the derivative is evaluated
        control_points : `TensorGrid` or `array-like`
            The points ``x_j`` controlling the deformation. They can
            be given either as a tensor grid or as a point array. In
            the latter case, its shape must be ``(N, n)``, where
            ``n`` is the dimension of the template space, and ``N``
            the number of ``alpha_j``, i.e. the size of (each
            component of) ``par_space``.
        discr_space : `DiscreteSpace`
            Space of the image grid of the template.
        kernel : `callable`
            Function to determine the kernel at the control points ``K(y_j)``
            The function must accept a real variable and return a real number.
        """

        super().__init__(alphas.space, control_points, discr_space, ft_kernel)
        self.discr_space = discr_space
        self.range_space = odl.ProductSpace(self.discr_space,
                                            self.discr_space.ndim)
        Operator.__init__(self, alphas.space, self.range_space, linear=True)
        self.alphas = alphas

    def _call(self, betas):
        """Implementation of ``self(betas)``.

        Parameters
        ----------
        betas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points. It should be in the same space as alpha.
        """

        return self.displacement_ft(betas)

    @property
    def adjoint(self):
        """Adjoint of the displacement derivative."""
        adj_op = DisplacementDerivativeAdjoint(
            self.alphas, self.control_points, self.discr_space, self.ft_kernel)

        return adj_op


class DisplacementDerivativeAdjoint(DisplacementDerivative):

    """Adjoint of the Displacement operator derivative.
    """
    def __init__(self, alphas, control_points, discr_space, ft_kernel):
        """Initialize a new instance.

        Parameters
        ----------
        alphas : `ProductSpaceVector`
            Displacement parameters in which the derivative is evaluated
        control_points : `TensorGrid` or `array-like`
            The points ``x_j`` controlling the deformation. They can
            be given either as a tensor grid or as a point array. In
            the latter case, its shape must be ``(N, n)``, where
            ``n`` is the dimension of the template space, and ``N``
            the number of ``alpha_j``, i.e. the size of (each
            component of) ``par_space``.
        discr_space : `DiscreteSpace`
            Space of the image grid of the template.
        kernel : `callable`
            Function to determine the kernel at the control points ``K(y_j)``
            The function must accept a real variable and return a real number.
        """

        super().__init__(alphas, control_points, discr_space, ft_kernel)

        # Switch domain and range
        self._domain, self._range = self._range, self._domain

    def _call(self, grad_func):
        """Implement ``self(func)```.

        Parameters
        ----------
        func : `DiscreteLpVector`
            Element of the image's gradient space.
        """
        # If control grid is not the image grid, the following result for
        # the ajoint is not right. Because the kernel matrix in fitting
        # term is not symmetric.
        return self.displacement_ft(grad_func)


class LinearizedDeformationOperator(Operator):

    """Operator mapping parameters to a fixed deformed template.

    This operator computes for a fixed template ``I`` the deformed
    template::

        v --> I(x + v(x))

    where the vector field ``v`` depends on the deformation parameters
    ``alpha_j`` as follows::

        v(x) = sum_j (K(x, x_j) * alpha_j)

    Here, ``K`` is the RKHS kernel matrix and each ``alpha_j`` is an
    element of ``R^n``.
    """

    def __init__(self, template):
        """Initialize a new instance.

        Parameters
        ----------
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        """
        # The operator maps from the parameter space to the template
        # (image) space.
        self.template = template
        self.domain_space = odl.ProductSpace(self.template.space,
                                             self.template.space.ndim)

        super().__init__(self.domain_space, self.template.space, linear=False)

    def _call(self, displacement):
        """Implementation of ``self(displacement)``.

        Parameters
        ----------
        displacement: `ProductSpaceVector`
            Linearized deformation parameters for image grid points.
        """
        image_pts = self.template.space.grid.points()
        image_pts += np.asarray(displacement).T

        return self.template.interpolation(image_pts.T, bounds_check=False)

    def linear_deform(self, template, displacement):
        """Implementation of ``self(template, displacement)``.

        Parameters
        ----------
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        displacement: `ProductSpaceVector`
            Linearized deformation parameters for image grid points.
        """
        image_pts = template.space.grid.points()
        image_pts += np.asarray(displacement).T
        return template.interpolation(image_pts.T, bounds_check=False)

    def derivative(self, displacement):
        """Frechet derivative of this operator in ``disp``.

        Parameters
        ----------
        displacement: `ProductSpaceVector`
            Linearized deformation parameters for image grid points.

        Returns
        -------
        deriv : `Operator`
        The derivative of this operator, evaluated at ``displacement``
        """
        deriv_op = LinearizedDeformationDerivative(self.template,
                                                   displacement)
        return deriv_op


class LinearizedDeformationDerivative(LinearizedDeformationOperator):

    """Frechet derivative of the Linearized deformation operator."""

    def __init__(self, template, displacement):
        """Initialize a new instance.

        Parameters
        ----------
        displacement: `ProductSpaceVector`
            Linearized deformation parameters for image grid points.
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        """

        super().__init__(template)

        self.template = template
        self.displacement = displacement

        Operator.__init__(self, self.displacement.space, self.template.space,
                          linear=False)

    def derivative(self, *args, **kwargs):
        raise NotImplementedError

    def _deform_grad(self, grad_f):
        """Compute the deformation of template's gradient.

        deformed template's gradient: \gradient I(x + v(x))

        Parameters
        ----------
        grad_f: `ProductSpaceVector`
            Gradient of the template, i.e. a vector field on the
            image domain
        """
        temp = [self.linear_deform(gf, self.displacement) for gf in grad_f]
        return self.displacement.space.element(temp)

    def _call(self, disp_like):
        """Implementation of ``self(disp_like)``.

        Parameters
        ----------
        disp_like: `ProductSpaceVector`
            Like the linearized deformation parameters for image grid points.
        """
        grad = odl.Gradient(self._range)
        template_grad = grad(self.template)
        def_grad = self._deform_grad(template_grad)

        result = def_grad[0] * disp_like[0].reshape(self.image_grid.shape,
                                                    order=self.range.order)
        for gf, d in zip(def_grad[1:], disp_like[1:]):
            result += gf * d.reshape(self.image_grid.shape,
                                     order=self.range.order)
        return result

    @property
    def adjoint(self):
        """Adjoint of the template deformation derivative."""
        adj_op = LinearizedDeformationDerivativeAdjoint(
            self.template, self.displacement)
        return adj_op


class LinearizedDeformationDerivativeAdjoint(LinearizedDeformationDerivative):

    """Adjoint of the template deformation operator derivative.

    TODO: write a bit more
    """
    def __init__(self, template, displacement):
        """Initialize a new instance.

        Parameters
        ----------
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        displacement: `ProductSpaceVector`
            Linearized deformation parameters for image grid points.
        """
        super().__init__(template, displacement)

        # Switch domain and range
        self._domain, self._range = self._range, self._domain

    def _call(self, func):
        """Implement ``self(func)```.

        Parameters
        ----------
        func : `DiscreteLpVector`
            Element of the image space
        """
        grad = odl.Gradient(self._domain)
        template_grad = grad(self.template)

        def_grad = self._deform_grad(template_grad)

        def_grad = template_grad.space.element([gf * func for gf in def_grad])

        return def_grad


class ShapeRegularizationFunctional(Operator):

    """Regularization functional for linear shape deformations.

    The shape regularization functional is given as

        S(alpha) = 1/2 * ||v(alpha)||^2 = 1/2 * alpha^T K alpha,

    where ``||.||`` is the norm in a reproducing kernel Hilbert space
    given by parameters ``alpha``. ``K`` is the kernel matrix.
    """
    # TODO: let user specify K

    def __init__(self, par_space, ft_kernel):
        """Initialize a new instance.

        Parameters
        ----------
        par_space : `ProductSpace`
            Parameter space of the deformations, i.e. the space of the
            ``alpha_k`` parametrizing the deformations
        kernel_op : `numpy.ndarray` or `Operator`
            The kernel matrix defining the norm. If an operator is
            given, it is called to evaluate the matrix-vector product
            of the kernel matrix with a given ``alpha``.
        """
        super().__init__(par_space, odl.RealNumbers(), linear=False)
#        if isinstance(kernel_op, Operator):
#            self._kernel_op = kernel_op
#        else:
#            self._kernel_op = odl.MatVecOperator(kernel_op)
        self.par_space = par_space
        self.ft_kernel = ft_kernel

    def _call(self, alphas):
        """Return ``self(alphas)``."""
        # TODO: how to improve

        # Compute the shape energy by fft
        ft_momenta = vectorial_ft_shape_op(alphas)
        stack = vectorial_ft_shape_op.inverse(self.ft_kernel * ft_momenta)
        return sum(s.inner(s.space.element(
                       np.asarray(a).reshape(-1, order=self.domain[0].order)))
                   for s, a in zip(stack, alphas)) / 2

#        # Compute the shape energy by matrix times vector
#        stack = [self._kernel_op(
#                     np.asarray(a).reshape(-1, order=self.domain[0].order))
#                 for a in alphas]
#        return sum(s.inner(s.space.element(
#                       np.asarray(a).reshape(-1, order=self.domain[0].order)))
#                   for s, a in zip(stack, alphas)) / 2

    def _gradient(self, alphas):
        """Return the gradient at ``alphas``.

        The gradient of the functional is given by

            grad(S)(alpha) = K alpha
        """
#        return self.domain.element([self._kernel_op(np.asarray(a).reshape(-1))
#                                    for a in alphas])
        pass

    def _gradient_ft(self, alphas):
        """Return the gradient at ``alphas``.

        The gradient of the functional is given by

            grad(S)(alpha) = K alpha.

        This is used for the n-D case: control grid = image grid.
        """
        ft_momenta = vectorial_ft_shape_op(alphas)
        return vectorial_ft_shape_op.inverse(self.ft_kernel * ft_momenta)
#        return (vectorial_ft_op_inverse(ft_displacement) /
#                self.par_space[0].cell_volume * 2.0 * np.pi)


class L2DataMatchingFunctional(Functional):

    """Basic data matching functional using the L2 norm.

    This functional computes::

        1/2 * ||f - g||_2^2

    for a given element ``g``.
    """

    def __init__(self, space, data):
        """Initialize a new instance.

        Parameters
        ----------
        space : `DiscreteLp` with exponent 2.0
            Space where the data is matched
        data : `DiscreteLp` element-like
            Data which is to be matched
        """
        if not (isinstance(space, odl.DiscreteLp) and space.exponent == 2.0):
            raise ValueError('not an L2 space.')
        super().__init__(space, linear=False)
        self.data = self.domain.element(data)

    def _call(self, x):
        """Return ``self(x)``."""
        return self.domain.dist(x, self.data)

    def gradient(self, x):
        """Return the gradient in the point ``x``."""
        return x - self.data

    def derivative(self, x):
        """Return the derivative in ``x``."""
        return self.gradient(x).T


# Kernel function for any dimensional
def gauss_kernel(x, sigma):
    return np.exp(-x ** 2 / (2 * sigma ** 2))


# Kernel function
def kernel(x):
    scaled = [xi ** 2 / (2 * sigma ** 2) for xi in x]
    return np.exp(-sum(scaled))


def gaussian_kernel_matrix(grid, sigma):
    """Return the kernel matrix for Gaussian kernel.

    The Gaussian kernel matrix ``K`` in ``n`` dimensions is defined as::

        k_ij = exp(- |x_i - x_j|^2 / (2 * sigma^2))

    where ``x_i, x_j`` runs through all grid points. The final matrix
    has size ``N x N``, where ``N`` is the total number of grid points.

    Parameters
    ----------
    grid : `RegularGrid`
        Grid where the control points are defined
    sigma : `float`
        Width of the Gaussian kernel
    """
    point_arrs = grid.points().T
    matrices = [parr[:, None] - parr[None, :] for parr in point_arrs]
    for mat in matrices:
        mat *= mat

    sq_sum = np.sqrt(np.sum(mat for mat in matrices))
    kernel_matrix = np.exp(-sq_sum / (2 * sigma ** 2))
    return kernel_matrix


def proj_noise(proj_data_shape, mu=0.0, sigma=0.1):
    """Produce white Gaussian noise for projections of n-D images.

       Produce white Gaussian noise for projections, with the same size
       as the number of projections.

       Parameters
       ----------
       proj_data_shape : shape of the projections
           Give the size of noise
       mu : Mean of the norm distribution
           The defalt is 0.0
       sigma : Standard deviation of the norm distribution
           The defalt is 0.1.
    """

    return np.random.normal(mu, sigma, proj_data_shape)


def SNR(signal, noise):
    """Compute the signal-to-noise ratio in dB.
    This compute::

    SNR = 10 * log10 (
        |signal - mean(signal)| / |noise - mean(noise)|)

    Parameters
    ----------
    signal : projection
    noise : white noise
    """
    if np.abs(np.asarray(noise)).sum() != 0:
        ave1 = np.sum(signal)/signal.size
        ave2 = np.sum(noise)/noise.size
        en1 = np.sqrt(np.sum((signal - ave1) * (signal - ave1)))
        en2 = np.sqrt(np.sum((noise - ave2) * (noise - ave2)))

        return 10.0 * np.log10(en1/en2)
    else:
        return 1000.

    return 10.0 * np.log10(en1/en2)


def padded_ft_op(space, padding_size):
    """Create zero-padding fft setting

    Parameters
    ----------
    space : the space needs to do FT
    padding_size : the percent for zero padding
    """
    padding_op = odl.ZeroPaddingOperator(
        space, [padding_size for _ in range(space.ndim)])
    shifts = [not s % 2 for s in space.shape]
    ft_op = odl.trafos.FourierTransform(
        padding_op.range, halfcomplex=False, shift=shifts)

    return ft_op * padding_op


def shape_kernel_ft(kernel):
    """Compute the n-D Fourier transform of the discrete kernel ``K``.

    Calculate the n-D Fourier transform of the discrete kernel ``K`` on the
    control grid points {y_i} to its reciprocal points {xi_i}.
    """

    # Create the array of kernel values on the grid points
    discretized_kernel = vspace.element(
        [cptsspace.element(kernel) for _ in range(cptsspace.ndim)])
    return vectorial_ft_shape_op(discretized_kernel)


def fitting_kernel_ft(kernel):
    """Compute the n-D Fourier transform of the discrete kernel ``K``.

    Calculate the n-D Fourier transform of the discrete kernel ``K`` on the
    image grid points {y_i} to its reciprocal points {xi_i}.

    """
    kspace = odl.ProductSpace(discr_space, discr_space.ndim)

    # Create the array of kernel values on the grid points
    discretized_kernel = kspace.element(
        [discr_space.element(kernel) for _ in range(discr_space.ndim)])
    return vectorial_ft_fit_op(discretized_kernel)


## Fix the sigma parameter in the kernel
#sigma = 5.0
#
## Discretization of the space, number of gridpoints for discretization
#m = 101
#
## Create 2-D discretization reconstruction space
## discr_space = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [m, m],
##                                dtype='float32', interp='linear')
#
## Create 3-D discretization reconstruction space
#discr_space = odl.uniform_discr(
#    [-0.5, -0.5, -0.5], [0.5, 0.5, 0.5], [m, m, m],
#    dtype='float32', interp='linear')
#
## Deformation space, number of gridpoints for deformation, usually n << m
#n = 101
#
## Create 2-D discretization space for control points
## cptsspace = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [n, n],
##                               dtype='float32', interp='linear')
#
## Create 3-D discretization space for control points
#cptsspace = odl.uniform_discr([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5], [n, n, n],
#                              dtype='float32', interp='linear')
#
## Create discretization space for vector field
#vspace = odl.ProductSpace(cptsspace, cptsspace.ndim)
#
## Create input function as Shepp-Logan phantom
#template = odl.util.shepp_logan(discr_space, modified=True)
#
## Create ground_truth function as Shepp-Logan phantom
#ground_truth = odl.util.shepp_logan(discr_space, modified=True)
#
## Create input function as disc phantom
## template = odl.util.disc_phantom(discr_space, smooth=True, taper=50.0)
## template.show('Template')
#
## Create ground_truth function as submarine phantom
## ground_truth = odl.util.submarine_phantom(
##    discr_space, smooth=True, taper=50.0)
## ground_truth.show('Ground Truth')
#
## Create 2-D projection domain
## detector_partition = odl.uniform_partition(-0.75, 0.75, 151)
#
## Create 3-D projection domain
#detector_partition = odl.uniform_partition(
#    [-0.75, -0.75], [0.75, 0.75], [151, 151])
#
## Create angle partition, Create projection directions
#angle_partition = odl.uniform_partition(0, np.pi/4, 2, nodes_on_bdry=True)
#
## Create 2-D parallel projection geometry
## geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)
#
## Create 3-D axis parallel projection geometry
#geometry = odl.tomo.Parallel3dAxisGeometry(angle_partition, detector_partition)

# Fix the sigma parameter in the kernel
sigma = 2.0

# Give the path of images
I0name = '../ddmatch/Example3 letters/c_highres.png'
I1name = '../ddmatch/Example3 letters/i_highres.png'

# Read the images to play with
I0 = plt.imread(I0name).astype('float')
I1 = plt.imread(I1name).astype('float')

# Do a downsampling
I0 = I0[::2, ::2]
I1 = I1[::2, ::2]

# Create 2-D discretization reconstruction space
# The size of the domain should be proportional to the given images
discr_space = odl.uniform_discr([-16, -16],
                                [16, 16], [128, 128],
                                dtype='float32', interp='linear')

# Create 2-D discretization space for control points
cptsspace = odl.uniform_discr([-16, -16], [16, 16], [128, 128],
                              dtype='float32', interp='linear')

# Create discretization space for vector field
vspace = odl.ProductSpace(cptsspace, cptsspace.ndim)

# Create the ground truth as the given image
ground_truth = discr_space.element(I0.T)

# Create the template as the given image
template = discr_space.element(I1.T)

# Give the number of directions
num_angles = 6

# Create the uniformly distributed directions
angle_partition = odl.uniform_partition(
    0, np.pi, num_angles, nodes_on_bdry=[(True, False)])

# Create 2-D projection domain
# The length should be 1.5 times of that of the reconstruction space
detector_partition = odl.uniform_partition(-24, 24, 192)

# Create 2-D parallel projection geometry
geometry = odl.tomo.Parallel2dGeometry(angle_partition,
                                       detector_partition)

# Create forward projection operator by X-ray transform
xray_trafo_op = odl.tomo.RayTransform(discr_space,
                                      geometry,
                                      impl='astra_cuda')

# Create projection data by given setting
proj_data = xray_trafo_op(ground_truth)

# Create white Gaussian noise
noise = 10.0 * proj_data.space.element(proj_noise(proj_data.shape))

# Compute the signal-to-noise ratio
snr = SNR(proj_data, noise)

# Output the signal-to-noise ratio
print('snr = {!r}'.format(snr))

# Create the noisy projection data
noise_proj_data = proj_data + noise

# Do the backprojection reconstruction
backproj = xray_trafo_op.adjoint(noise_proj_data)

# FFT setting for regularization shape term, 1 means 100% padding
padding_size = 1.0
padded_ft_shape_op = padded_ft_op(cptsspace, padding_size)
vectorial_ft_shape_op = odl.DiagonalOperator(
    *([padded_ft_shape_op] * cptsspace.ndim))

# FFT setting for data matching term, 1 means 100% padding
padded_ft_fit_op = padded_ft_op(discr_space, padding_size)
vectorial_ft_fit_op = odl.DiagonalOperator(
    *([padded_ft_fit_op] * cptsspace.ndim))

# Initialize deformation vector field
momenta = vspace.zero()

# Compute Fourier trasform of the kernel function in data matching term
ft_kernel_fitting = fitting_kernel_ft(kernel)

# Compute Fourier trasform of the kernel function in shape regularization term
ft_kernel_shape = shape_kernel_ft(kernel)

# Create displacement operator
displacement_op = DisplacementOperator(vspace, cptsspace.grid,
                                       discr_space, ft_kernel_fitting)

# Compute the displacement at momenta
displ = displacement_op(momenta)

# Create linearized deformation operator
linear_deform_op = LinearizedDeformationOperator(template)

# Compute the deformed template
deformed_template = linear_deform_op(displ)

# Create X-ray transform operator
proj_deformed_template = xray_trafo_op(deformed_template)

# Create L2 data matching (fitting) term
l2_data_fit_func = L2DataMatchingFunctional(xray_trafo_op.range,
                                            noise_proj_data)

# Composition of the L2 fitting term with three operators
data_fitting_term = l2_data_fit_func * xray_trafo_op * linear_deform_op * displacement_op

# Compute the kernel matrix for the method without Fourier transform
# If the dimension is too large, it could cause MemoryError
# kernelmatrix = gaussian_kernel_matrix(cptsspace.grid, sigma)

# Compute the gradient of shape regularization term
shape_func = ShapeRegularizationFunctional(vspace, ft_kernel_shape)

# Shape regularization parameter, should be nonnegtive
lambda_shape = 0.0001

# Step size for the gradient descent method
eta = 0.05

# Maximum iteration number
niter = 500

# Test time, set starting time
start = time.clock()

# Create the memory for energy in each iteration
E = []
kE = len(E)
E = np.hstack((E, np.zeros(niter)))

# Iterations for updating alphas
for i in range(niter):

    # Compute the gradient for the shape term by Fourier transform
    grad_shape_func = shape_func._gradient_ft(momenta)

    E[i+kE] = 0.5 * lambda_shape * sum(s.inner(s.space.element(np.asarray(a).reshape(-1, order=vspace[0].order)))
        for s, a in zip(grad_shape_func, momenta)) + data_fitting_term(momenta)

    # Compute the gradient for data fitting term
    grad_data_fitting_term = data_fitting_term.gradient(momenta)

    # Update momenta
    momenta -= eta * (
        lambda_shape * grad_shape_func + grad_data_fitting_term)

    # Show the middle reconstrcted results
    if (i+1) % 100 == 0:
        print(i+1)
#        displ = displacement_op(momenta)
#        deformed_template = linear_deform_op(displ)
#        deformed_template.show(
#            title='Reconstruction Image, iters = {!r}, eta 200'.format(i+1))

# Test time, set end time
end = time.clock()

# Output the computational time
print(end - start)

# Compute the projections of the reconstructed image
displ = displacement_op(momenta)
deformed_template = linear_deform_op(displ)
rec_proj_data = xray_trafo_op(deformed_template)

# Plot the results of interest
plt.figure(1, figsize=(20, 10))
plt.clf()

plt.subplot(2, 3, 1)
plt.imshow(I0, cmap='bone', vmin=I0.min(), vmax=I0.max())
plt.colorbar()
plt.title('Ground truth')

plt.subplot(2, 3, 2)
plt.imshow(I1, cmap='bone', vmin=I1.min(), vmax=I1.max())
plt.colorbar()
plt.title('Template')

plt.subplot(2, 3, 3)
# plt.imshow(np.asarray(backproj).T, cmap='bone',
#            vmin=np.asarray(backproj).min(), vmax=np.asarray(backproj).max())
# plt.colorbar()
# plt.title('Backprojection')
plt.title('stepsize = {!r}, $\sigma$ = {!r}'.format(eta, sigma))
plt.plot(E)
plt.xlabel('Iteration number')
plt.ylabel('Energy')
# plt.gca().axes.yaxis.set_ticklabels(['0']+['']*8)
plt.gca().axes.yaxis.set_ticklabels([])
plt.grid(True)

plt.subplot(2, 3, 4)
plt.imshow(np.asarray(deformed_template).T, cmap='bone',
           vmin=np.asarray(deformed_template).min(),
           vmax=np.asarray(deformed_template).max())
plt.colorbar()
plt.title('Reconstructed image by {!r} iters, {!r} projs'.format(niter, num_angles))

plt.subplot(2, 3, 5)
plt.plot(np.asarray(proj_data)[0], 'b', np.asarray(noise_proj_data)[0], 'r')
plt.title('Theta=0, blue: truth_data, red: noisy_data, SNR = {:.3}dB'.format(snr))
plt.gca().axes.yaxis.set_ticklabels([])
plt.axis([0, 191, -17, 32])

plt.subplot(2, 3, 6)
plt.plot(np.asarray(proj_data)[0], 'b', np.asarray(rec_proj_data)[0], 'r')
plt.title('Theta=0, blue: truth_data, red: rec result')
plt.gca().axes.yaxis.set_ticklabels([])
plt.axis([0, 191, -17, 32])
