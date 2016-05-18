# -*- coding: utf-8 -*-
"""
Example of creating an operator that acts as a linear transformation.
"""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super
from numbers import Number
from odl.operator.operator import Operator, OperatorComp
import odl
import numpy as np
from functools import partial
from time import clock
from scipy import signal


# Computing deformation at point x
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


class TemplateDeformationOperator(Operator):

    """Operator mapping parameters to a fixed deformed template.

    This operator computes for a fixed template ``I`` the deformed
    template::

        I --> I(x + v(x))

    where the vector field ``v`` depends on the deformation parameters
    ``alpha_j`` as follows::

        v(x) = sum_j (K(x, x_j) * alpha_j)

    Here, ``K`` is the RKHS kernel matrix and each ``alpha_j`` is an
    element of ``R^n``.
    """

    def __init__(self, par_space, control_points, template, kernel,
                 kernel_func=None, **kwargs):
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
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        kernel_func : `callable`
            Function to determine the kernel matrix ``K(x, x_j)``
            as ``kernel_func(|x - x_j|) * eye(n)``. The function must
            accept a real variable and return a real number.
        cache_kernel_matrix : `bool`
            If `True`, store the kernel matrix after it is first
            computed.
            Default: `False`
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

        if not isinstance(template, odl.DiscreteLpVector):
            raise TypeError('expected DiscreteLpVector as template, got {!r}.'
                            ''.format(template))

        if par_space.size != template.space.ndim:
            raise ValueError('dimensions of product space and template space '
                             'do not match ({} != {})'
                             ''.format(par_space.size, template.space.ndim))

        # The operator maps from the parameter space to the template
        # (image) space.
        super().__init__(par_space, template.space, linear=False)
        self.template = template

        if not callable(kernel_func):
            raise TypeError('kernel function {!r} is not callable.'
                            ''.format(kernel_func))

        self.kernel_func = kernel_func

        self.kernel = kernel

        if not isinstance(control_points, odl.TensorGrid):
            self._control_pts = np.asarray(control_points)
            if self._control_pts.shape != (self.num_contr_pts, self.ndim):
                raise ValueError(
                    'expected control point array of shape {}, got {}.'
                    ''.format((self.num_contr_pts, self.ndim),
                              self.control_points.shape))
        else:
            self._control_pts = control_points

        self._cache_kernel_matrix = bool(kwargs.pop('cache_kernel_matrix',
                                                    False))
        self._kernel_matrix = None
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
        return self.template.space.grid

    def kernel_matrix(self, c_idx=None):
        """Compute the kernel matrix ``K``.

        Calculate the matrix of differences |x_i - y_j| for the selected
        control points x_i and image grid points y_j.

        Parameters
        ----------
        c_idx : index expression
            Choose control points (= rows) according to the given
            indices. If the control points are given as an array,
            one-dimensional indices must be used. Otherwise, an
            ``n``-dimensional index expression is needed, where
            ``n`` is ``self.ndim``.
        """
        if self._kernel_matrix is not None:
            return self._kernel_matrix

        if c_idx is None:
            if self.contr_pts_is_grid:
                contr_pt_arrs = self.control_points.points().T
            else:
                contr_pt_arrs = self.control_points.T
        else:
            # Select control point according to c_idx
            contr_pt_arrs = self.control_points[c_idx]
            if isinstance(contr_pt_arrs, odl.TensorGrid):
                contr_pt_arrs = contr_pt_arrs.points().T
            else:
                contr_pt_arrs = np.atleast_2d(contr_pt_arrs.T)

        image_pt_arrs = self.image_grid.points().T

        # Calculate the matrix of differences |x_i - y_j| for the selected
        # control points x_i and image grid points y_j
        matrices = [px[:, None] - py[None, :]
                    for px, py in zip(contr_pt_arrs, image_pt_arrs)]
        for mat in matrices:
            mat *= mat

        abs_diff = sum(mat for mat in matrices)
        np.sqrt(abs_diff, out=abs_diff)
        kernel_matrix = self.kernel_func(abs_diff)

        if self._cache_kernel_matrix:
            self._kernel_matrix = kernel_matrix

        return kernel_matrix

    def kernel_2dfft(self):
        """Compute the 2D Fourier transform of the discrete kernel ``K``.

        Calculate the 2D Fourier transform of the discrete kernel ``K`` on the
        grid points {x_i} to its reciprocal points {xi_i}.

        """
        kspace = odl.ProductSpace(self.template.space, 2)

        # Create the array of kernel values on the grid points
        discretized_kernel1 = self.template.space.element(self.kernel)
        discretized_kernel2 = self.template.space.element(self.kernel)
        discretized_kernel = kspace.element([discretized_kernel1,
                                             discretized_kernel2])

        # Make the Fourier transform operator on this space. The range is
        # calculated automatically. The default backend is numpy.fft. Also
        # ensure the origin is in the reciprocal grid.
        shifts = [not s % 2 for s in self.template.space.shape]
        ft_op = odl.trafos.FourierTransform(
            self.template.space, halfcomplex=False, shift=shifts)
        vectorial_ft_op = odl.ProductSpaceOperator([[ft_op, 0],
                                                    [0, ft_op]])

        # Step 1
        # Calculate the Fourier transform of the vectorial kernel (step 1)
        ft_kernel = vectorial_ft_op(discretized_kernel)

        return ft_kernel

    def kernel_2dfft_zero_padding(self):
        """Compute the 2D Fourier transform of the discrete kernel ``K``.

        Calculate the 2D Fourier transform of the discrete kernel ``K`` on the
        grid points {x_i} to its reciprocal points {xi_i}.

        """
        kspace = odl.ProductSpace(self.template.space, 2)

        # Create the array of kernel values on the grid points
        discretized_kernel1 = self.template.space.element(self.kernel)
        discretized_kernel2 = self.template.space.element(self.kernel)
        discretized_kernel = kspace.element([discretized_kernel1,
                                             discretized_kernel2])

        padding_op = odl.ZeroPaddingOperator(self.template.space, [1, 1])
        shifts = [not s % 2 for s in self.template.space.shape]

        ft_op = odl.trafos.FourierTransform(
            padding_op.range, halfcomplex=False, shift=shifts)

        padded_ft_op = ft_op * padding_op
        vectorial_ft_op = odl.ProductSpaceOperator([[padded_ft_op, 0],
                                                    [0, padded_ft_op]])

        ft_kernel = vectorial_ft_op(discretized_kernel)

        return ft_kernel

    def displacement(self, alphas):
        """Calculate the translation at point x.

        displacement: v(x)

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points. Note that here N = M is not neccessary.
        """
        kernel_mat = self.kernel_matrix()

        print(kernel_mat)
        print(sum(kernel_mat))

        displacement = [kernel_mat.T.dot(np.asarray(a).reshape(
            -1, order=alphas.space[0].order))
                for a in alphas]

        return displacement

    def displacement_2dfft(self, alphas):
        """Calculate the translation at point x by 2D FFT.

        displacement: v(x)

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points. Note that here N = M.
        kernel: Gaussian kernel fuction
        """
        shifts = [not s % 2 for s in self.template.space.shape]
        ft_op = odl.trafos.FourierTransform(
            self.template.space, halfcomplex=False, shift=shifts)
        vectorial_ft_op = odl.ProductSpaceOperator([[ft_op, 0],
                                                    [0, ft_op]])

        vectorial_ft_op_inverse = odl.ProductSpaceOperator(
            [[ft_op.inverse, 0],
             [0, ft_op.inverse]])

        ft_momenta = vectorial_ft_op(alphas)
        ft_displacement = self.kernel_2dfft(self.kernel)*ft_momenta
        return vectorial_ft_op_inverse(ft_displacement)

    def displacement_2dfft_zero_padding(self, alphas):
        """Calculate the translation at point x by 2D FFT.

        displacement: v(x)

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points. Note that here N = M.
        kernel: Gaussian kernel fuction
        """
        padding_op = odl.ZeroPaddingOperator(self.template.space, [1, 1])
        shifts = [not s % 2 for s in self.template.space.shape]

        ft_op = odl.trafos.FourierTransform(
            padding_op.range, halfcomplex=False, shift=shifts)

        padded_ft_op = ft_op * padding_op
        vectorial_ft_op = odl.ProductSpaceOperator([[padded_ft_op, 0],
                                                    [0, padded_ft_op]])

        vectorial_ft_op_inverse = odl.ProductSpaceOperator(
            [[padded_ft_op.inverse, 0],
             [0, padded_ft_op.inverse]])
        ft_momenta = vectorial_ft_op(alphas)
        ft_displacement = self.kernel_2dfft_zero_padding()*ft_momenta
        return vectorial_ft_op_inverse(ft_displacement)

    def displacement_2dconv(self, alphas):
        """Calculate the translation at point x by 2D FFT.

        displacement: v(x)

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points. Note that here N = M.
        kernel: Gaussian kernel fuction
        """

        kspace = odl.ProductSpace(self.template.space, 2)

        # Create the array of kernel values on the grid points
        discretized_kernel1 = self.template.space.element(self.kernel)
        discretized_kernel2 = self.template.space.element(self.kernel)

        displacement1 = signal.convolve2d(alphas[0], discretized_kernel1,
                                          mode='same', boundary='symm',
                                          fillvalue=0)
        displacement2 = signal.convolve2d(alphas[1], discretized_kernel2,
                                          mode='same', boundary='symm',
                                          fillvalue=0)

        return kspace.element([displacement1, displacement2])

    def deform(self, template, alphas):
        """Compute the deformation of the template.

        deformed template: I(x + v(x))

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.
        """
        image_pts = self.image_grid.points()
        displacement = self.displacement(alphas)
        image_pts += np.asarray(displacement).T

        return template.interpolation(image_pts.T, bounds_check=False)

    def deform_2dfft(self, template, alphas):
        """Compute the deformation of the template by Fourier transform.

        deformed template: I(x + v(x))

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.
        """
        image_pts = template.space.grid.points()
        displacement = self.displacement_2dfft(alphas, self.kernel)

        # scaling
        displacement = displacement * 2.0 * np.pi / template.cell_volume

        print('displacement 2dfft')
        print(displacement)

        displacement.show('Displacement_2dfft')

        image_pts += np.asarray(displacement).T

        return template.interpolation(image_pts.T, bounds_check=False)

    def deform_2dconv(self, template, alphas):
        """Compute the deformation of the template by Fourier transform.

        deformed template: I(x + v(x))

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.
        """
        image_pts = template.space.grid.points()
        displacement = self.displacement_2dconv(alphas, self.kernel)

        print('displacement 2dconv')
        print(displacement)

        displacement.show('Displacement_2dconv')

        image_pts += np.asarray(displacement).T

        return template.interpolation(image_pts.T, bounds_check=False)

    def deform_2dfft_zero_padding(self, template, alphas):
        """Compute the deformation of the template by Fourier transform.

        deformed template: I(x + v(x))

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.
        """
        image_pts = template.space.grid.points()
        displacement = self.displacement_2dfft_zero_padding(alphas)
        # Rescale with the 1/cell_volume
        displacement /= discr_space.cell_volume
        displacement *= 2 * np.pi

#        print('displacement 2dfft_zero_padding')
#        print(displacement)

#        displacement.show('Displacement_2dfft_zero_padding')

        image_pts += np.asarray(displacement).T

#        return template.space.element(template.interpolation(
#                                            image_pts.T, bounds_check=False))
        return template.interpolation(image_pts.T, bounds_check=False)

    def _call(self, alphas, out):
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
        # TODO: Lazy hack, optimize later
        # out[:] = self.deform(self.template, alphas)# mat-vec method
        out[:] = self.deform_2dfft_zero_padding(self.template, alphas)

    def derivative(self, alphas, **kwargs):
        """Frechet derivative of this operator in ``alphas``.

        Parameters
        ----------
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.
        cache_deformed_grad : `bool`
            If `True`, the deformed gradient of the template is
            cached.
            Default: `True`

        Returns
        -------
        deriv : `Operator`
            The derivative of this operator, evaluated at ``alphas``
        """
        cache_deformed_grad = kwargs.pop('cache_deformed_grad', True)
        deriv_op = TemplateDeformationDerivative(
            alphas, self.control_points, self.template, self.kernel,
            self.kernel_func, cache_kernel_matrix=self._cache_kernel_matrix,
            cache_deformed_grad=cache_deformed_grad)

        # Share the cached stuff
        deriv_op._kernel_matrix = self._kernel_matrix
        return deriv_op


class TemplateDeformationDerivative(TemplateDeformationOperator):

    """Frechet derivative of the template deformation operator."""

    def __init__(self, alphas, control_points, template, kernel,
                 kernel_func=None, **kwargs):
        """Initialize a new instance.

        Parameters
        ----------
        alphas : `ProductSpaceVector`
            Deformation parameters in which the derivative is evaluated
        control_points : `TensorGrid` or `array-like`
            The points ``x_j`` controlling the deformation. They can
            be given either as a tensor grid or as a point array. In
            the latter case, its shape must be ``(N, n)``, where
            ``n`` is the dimension of the template space, and ``N``
            the number of ``alpha_j``, i.e. the size of (each
            component of) ``par_space``.
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        kernel_func : `callable`
            Function to determine the kernel matrix ``K(x, x_j)``
            as ``kernel_func(|x - x_j|) * eye(n)``. The function must
            accept a real variable and return a real number.
        cache_kernel_matrix : `bool`
            If `True`, store the kernel matrix after it is first
            computed.
            Default: `False`
        cache_deformed_grad : `bool`
            If `True`, the deformed gradient of the template is
            cached.
            Default: `True`
        """
        super().__init__(alphas.space, control_points, template, kernel,
                         kernel_func, **kwargs)
        Operator.__init__(self, alphas.space, template.space, linear=True)
        self.alphas = alphas
        cache_deformed_grad = kwargs.pop('cache_deformed_grad', True)
        self._cache_deformed_grad = bool(cache_deformed_grad)
        self._deformed_grad = None

    def derivative(self, *args, **kwargs):
        raise NotImplementedError

    def _deform_grad(self, grad_f, alphas, out=None):
        """Compute the deformation of template's gradient.

        deformed template's gradient: \gradient I(x + v(x))

        Parameters
        ----------
        grad_f: `ProductSpaceVector`
            Gradient of the template, i.e. a vector field on the
            image domain
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.
        out : `ProductSpaceVector`, optional
            Element where the result is stored

        Returns
        -------
        out : `ProductSpaceVector`
            The deformed template gradient. If ``out`` was given, the
            returned object is a reference to it.
        """
        grad_space = odl.ProductSpace(self.range, self.ndim)
        if out is None:
            out = grad_space.element()

        for gf, oi in zip(grad_f, out):
            # TODO: do deformation in place
            # oi[:] = self.deform(gf, alphas)
            oi[:] = self.deform_2dfft_zero_padding(gf, alphas)

        return out

    def _call(self, betas, out):
        """Implementation of ``self(betas, out)``.

        Parameters
        ----------
        betas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.

        out : `DiscreteLpVector`
            Element where the result is stored
        """
        # Use cached value if possible
        if self._deformed_grad is not None:
            def_grad = self._deformed_grad
        else:
            grad = odl.Gradient(self.range)
            template_grad = grad(self.template)
            def_grad = self._deform_grad(template_grad, self.alphas,
                                         out=template_grad)
            if self._cache_deformed_grad:
                self._deformed_grad = def_grad

        displ = self.displacement(betas)

        out[:] = def_grad[0] * displ[0].reshape(self.image_grid.shape,
                                                order=self.range.order)
        for gf, d in zip(def_grad[1:], displ[1:]):
            out += gf * d.reshape(self.image_grid.shape,
                                  order=self.range.order)

    @property
    def adjoint(self):
        """Adjoint of the template deformation derivative."""
        adj_op = TemplateDeformationDerivativeAdjoint(
            self.alphas, self.control_points, self.template, self.kernel,
            self.kernel_func, cache_kernel_matrix=self._cache_kernel_matrix,
            cache_deformed_grad=self._cache_deformed_grad)

        # Share the cached stuff
        adj_op._kernel_matrix = self._kernel_matrix
        adj_op._deformed_grad = self._deformed_grad
        return adj_op


class TemplateDeformationDerivativeAdjoint(TemplateDeformationDerivative):

    """Adjoint of the template deformation operator derivative.

    TODO: write a bit more
    """
    def __init__(self, alphas, control_points, template, kernel,
                 kernel_func=None, **kwargs):
        """Initialize a new instance.

        Parameters
        ----------
        alphas : `ProductSpaceVector`
            Deformation parameters in which the derivative is evaluated
        control_points : `TensorGrid` or `array-like`
            The points ``x_j`` controlling the deformation. They can
            be given either as a tensor grid or as a point array. In
            the latter case, its shape must be ``(N, n)``, where
            ``n`` is the dimension of the template space, and ``N``
            the number of ``alpha_j``, i.e. the size of (each
            component of) ``par_space``.
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        kernel_func : `callable`
            Function to determine the kernel matrix ``K(x, x_j)``
            as ``kernel_func(|x - x_j|) * eye(n)``. The function must
            accept a real variable and return a real number.
        cache_kernel_matrix : `bool`
            If `True`, store the kernel matrix after it is first
            computed.
            Default: `False`
        cache_deformed_grad : `bool`
            If `True`, the deformed gradient of the template is
            cached.
            Default: `True`
        """
        super().__init__(alphas, control_points, template, kernel, kernel_func,
                         **kwargs)

        # Switch domain and range
        self._domain, self._range = self._range, self._domain

    def _call(self, func, out):
        """Implement ``self(func, out)```.

        Parameters
        ----------
        func : `DiscreteLpVector`
            Element of the image space
        out : `ProductSpaceVector`
            Element of the parameter space where the result is
            written into
        """
        # Use cached value if possible
        if self._deformed_grad is not None:
            def_grad = self._deformed_grad
        else:
            grad = odl.Gradient(self.domain)
            template_grad = grad(self.template)
            def_grad = self._deform_grad(template_grad, self.alphas,
                                         out=template_grad)
            if self._cache_deformed_grad:
                self._deformed_grad = def_grad

        for gf in def_grad:
            gf = gf * func
        out = self.displacement_2dfft_zero_padding(def_grad)
        out *= self.domain.cell_volume


#        for a, gf in zip(out, def_grad):
#            self.kernel_matrix().dot(
#                b=np.asarray(gf * func).reshape(-1, order=self.domain.order),
#                out=np.asarray(a).reshape(-1, order=self.range[0].order))
#        out *= self.domain.cell_volume


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

#    def kernel_2dfft_zero_padding(self):
#        """Compute the 2D Fourier transform of the discrete kernel ``K``.
#
#        Calculate the 2D Fourier transform of the discrete kernel ``K`` on the
#        grid points {y_i} to its reciprocal points {xi_i}.
#
#        """
#        kspace = odl.ProductSpace(discr_space, 2)
#
#        # Create the array of kernel values on the grid points
#        discretized_kernel1 = discr_space.element(self.kernel)
#        discretized_kernel2 = discr_space.element(self.kernel)
#        discretized_kernel = kspace.element([discretized_kernel1,
#                                             discretized_kernel2])
#
#        padding_op = odl.ZeroPaddingOperator(discr_space, [1, 1])
#        shifts = [not s % 2 for s in discr_space.shape]
#
#        ft_op = odl.trafos.FourierTransform(
#            padding_op.range, halfcomplex=False, shift=shifts)
#
#        padded_ft_op = ft_op * padding_op
#        vectorial_ft_op = odl.ProductSpaceOperator([[padded_ft_op, 0],
#                                                    [0, padded_ft_op]])
#        ft_kernel = vectorial_ft_op(discretized_kernel)
#        return ft_kernel

    def displacement_2dfft_zero_padding(self, alphas):
        """Calculate the inverse translation at point y by 2D FFT.

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
        padding_op = odl.ZeroPaddingOperator(discr_space, [1, 1])
        shifts = [not s % 2 for s in discr_space.shape]

        ft_op = odl.trafos.FourierTransform(
            padding_op.range, halfcomplex=False, shift=shifts)

        padded_ft_op = ft_op * padding_op
        vectorial_ft_op = odl.ProductSpaceOperator([[padded_ft_op, 0],
                                                    [0, padded_ft_op]])

        vectorial_ft_op_inverse = odl.ProductSpaceOperator(
            [[padded_ft_op.inverse, 0],
             [0, padded_ft_op.inverse]])

        ft_momenta = vectorial_ft_op(alphas)
        ft_displacement = self.ft_kernel * ft_momenta
        return vectorial_ft_op_inverse(ft_displacement)
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

        return self.displacement_2dfft_zero_padding(alphas)

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

        return self.displacement_2dfft_zero_padding(betas)

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
        # TODO: How to deal with the above problem for general case?
        return self.displacement_2dfft_zero_padding(grad_func)


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

    def __init__(self, par_space, kernel_op):
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
        if isinstance(kernel_op, Operator):
            self._kernel_op = kernel_op
        else:
            self._kernel_op = odl.MatVecOperator(kernel_op)
        self.par_space = par_space

    def _call(self, alphas):
        """Return ``self(alphas)``."""
        # TODO: add out parameter
        stack = [self._kernel_op(
                     np.asarray(a).reshape(-1, order=self.domain[0].order))
                 for a in alphas]
        return sum(s.inner(s.space.element(
                       np.asarray(a).reshape(-1, order=self.domain[0].order)))
                   for s, a in zip(stack, alphas)) / 2

    def _gradient(self, alphas):
        """Return the gradient at ``alphas``.

        The gradient of the functional is given by

            grad(S)(alpha) = K alpha
        """
        return self.domain.element([self._kernel_op(np.asarray(a).reshape(-1))
                                    for a in alphas])

    def _gradient_2dfft_zero_padding(self, ft_momenta, ft_kernel):
        """Return the gradient at ``alphas``.

        The gradient of the functional is given by

            grad(S)(alpha) = K alpha.

        This is used for the 2D case: control grid = image grid.
        """
        ft_displacement = ft_kernel * ft_momenta
        return ft_displacement
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


class FittingTermFunctional(Functional):

    """Fitting term functional for shape-based image reconstruction.

    The fitting term functional is given as

        T(alpha) = 1/2 * ||T(I(. + \sum_1^N K(., x_j)alpha_j)) - g||_2^2,

    where ``||.||_2`` is the L^2 norm given by parameters ``alpha``. The T,
    in the right hand side, denotes Radon transform. ``g`` is the
    mearsured data and ``K`` is the kernel matrix.
    """

    def __init__(self, par_space, kernel_op, data):
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
        self.data = data
        super().__init__(par_space, odl.RealNumbers(), linear=False)
        if isinstance(kernel_op, Operator):
            self._kernel_op = kernel_op
        else:
            self._kernel_op = odl.MatVecOperator(kernel_op)

    def _call(self, alphas):
        """Return ``self(alphas)``."""
        stack = [self._kernel_op(np.asarray(a).reshape(-1)) for a in alphas]
        return sum(s.inner(s.space.element(np.asarray(a).reshape(-1)))
                   for s, a in zip(stack, alphas)) / 2

    # Computing the gradient for the fitting term
    def grad_fitting_term(self, deformed_grad, backproj_diff, kernel_matrix):

        deformed_grad = [gf * backproj_diff for gf in deformed_grad]
        return [py.dot(px.T) for px, py in zip(deformed_grad, kernel_matrix)]

    # Computing the gradient for the fitting term using zero-padding 2dfft

    def grad_fitting_term_2dfft_zero_padding(self, template, deformed_grad,
                                             backproj_diff, kernel):
        """Return the gradient for the fitting term, computing with zero-padding 2dfft.

        The gradient of the functional is given by

            grad(S)(alpha) = K * (T*^(T(deformed_template) - g) *
            grad_deformed_template)
        """
        deformed_grad = [gf * backproj_diff for gf in deformed_grad]

        padding_op = odl.ZeroPaddingOperator(template.space, [1, 1])
        shifts = [not s % 2 for s in template.space.shape]

        ft_op = odl.trafos.FourierTransform(
            padding_op.range, halfcomplex=False, shift=shifts)

        padded_ft_op = ft_op * padding_op
        vectorial_ft_op = odl.ProductSpaceOperator([[padded_ft_op, 0],
                                                    [0, padded_ft_op]])

        vectorial_ft_op_inverse = odl.ProductSpaceOperator(
            [[padded_ft_op.inverse, 0],
             [0, padded_ft_op.inverse]])

        ft_momenta = vectorial_ft_op(deformed_grad)
        kspace = odl.ProductSpace(template.space, 2)

        # Create the array of kernel values on the grid points
        discretized_kernel1 = template.space.element(kernel)
        discretized_kernel2 = template.space.element(kernel)
        discretized_kernel = kspace.element([discretized_kernel1,
                                             discretized_kernel2])

        ft_kernel = vectorial_ft_op(discretized_kernel)
        ft_displacement = ft_kernel * ft_momenta
        return vectorial_ft_op_inverse(ft_displacement)


# Fix the sigma parameter in the kernel
sigma = 0.2


# Kernel function for any dimensional
def gauss_kernel(x, sigma):
    return np.exp(-x ** 2 / (2 * sigma ** 2))


# Kernel function
def kernel(x):
    scaled = [xi ** 2 / (2 * sigma ** 2) for xi in x]
    return np.exp(-sum(scaled))

kernel1 = partial(gauss_kernel, sigma=sigma)


# Produce noise for projections of 2D images
def proj_noise(x, y, mu=0.0, sigma=1.0):
    return sigma * np.random.rand(x, y) + mu


# Compute Signal-to-Noise Ratio
def SNR(signal, noise):
    ave1 = np.sum(signal)/signal.size
    ave2 = np.sum(noise)/noise.size
    en1 = np.sqrt(np.sum((signal - ave1) * (signal - ave1)))
    en2 = np.sqrt(np.sum((noise - ave2) * (noise - ave2)))
    return 10.0 * np.log10(en1/en2)

# Discretization of the space
m = 101  # Number of gridpoints for discretization
discr_space = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [m, m],
                                dtype='float32', interp='linear')

# Deformation space
n = 101  # Number of gridpoints for deformation, usually n << m
cptssapce = odl.uniform_discr([-0.5, -0.5], [0.5, 0.5], [n, n],
                              dtype='float32', interp='linear')
vspace = odl.ProductSpace(cptssapce, 2)

# Create input function as disc phantom
template = odl.util.disc_phantom(discr_space, smooth=True, taper=50.0)
template.show('Original template')

# Create target function as submarine phantom
target = odl.util.submarine_phantom(discr_space, smooth=True, taper=50.0)
target.show('Target')

# Create projection domain
detector_partition = odl.uniform_partition(-0.75, 0.75, 151)

# Create projection directions
angle_interval = odl.Interval(0, np.pi)
angle_grid = odl.TensorGrid([0, np.pi/6, np.pi/3, np.pi/2, np.pi*2/3, np.pi*5/6])
angle_partition = odl.RectPartition(angle_interval, angle_grid)

# Create 2D parallel projection geometry
geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)

# Create forward projections by Radon transform
xray_trafo_op = odl.tomo.RayTransform(discr_space, geometry, impl='astra_cuda')
proj_data = xray_trafo_op(target)

# Create white Gaussian noise
noise = 0.1 * proj_data.space.element(proj_noise(proj_data.shape[0],
                                                 proj_data.shape[1]))

# Compute Signal-to-Noise Ratio
print(SNR(proj_data, noise))

# Create noisy projections, noise ~ (0, 0.1)
noise_proj_data = proj_data + noise

# proj_data_template = xray_trafo_op(template)
backproj = xray_trafo_op.adjoint(noise_proj_data)
backproj.show('backprojection')

# FFT setting for shape
padding_op = odl.ZeroPaddingOperator(vspace[0], [1, 1])
shifts = [not s % 2 for s in vspace[0].shape]

ft_op = odl.trafos.FourierTransform(
    padding_op.range, halfcomplex=False, shift=shifts)

padded_ft_op = ft_op * padding_op
vectorial_ft_op = odl.ProductSpaceOperator([[padded_ft_op, 0],
                                            [0, padded_ft_op]])

vectorial_ft_op_inverse = odl.ProductSpaceOperator([[padded_ft_op.inverse, 0],
                                                    [0, padded_ft_op.inverse]])


def shape_kernel_2dfft_zero_padding(kernel):
    """Compute the 2D Fourier transform of the discrete kernel ``K``.

    Calculate the 2D Fourier transform of the discrete kernel ``K`` on the
    control grid points {y_i} to its reciprocal points {xi_i}.

    """

    # Create the array of kernel values on the grid points
    discretized_kernel1 = vspace[0].element(kernel)
    discretized_kernel2 = vspace[0].element(kernel)
    discretized_kernel = vspace.element([discretized_kernel1,
                                         discretized_kernel2])

    ft_kernel = vectorial_ft_op(discretized_kernel)
    return ft_kernel


def fitting_kernel_2dfft_zero_padding(kernel):
    """Compute the 2D Fourier transform of the discrete kernel ``K``.

    Calculate the 2D Fourier transform of the discrete kernel ``K`` on the
    image grid points {y_i} to its reciprocal points {xi_i}.

    """
    kspace = odl.ProductSpace(discr_space, 2)

    # Create the array of kernel values on the grid points
    discretized_kernel1 = discr_space.element(kernel)
    discretized_kernel2 = discr_space.element(kernel)
    discretized_kernel = kspace.element([discretized_kernel1,
                                         discretized_kernel2])

    padding_op = odl.ZeroPaddingOperator(discr_space, [1, 1])
    shifts = [not s % 2 for s in discr_space.shape]

    ft_op = odl.trafos.FourierTransform(
        padding_op.range, halfcomplex=False, shift=shifts)

    padded_ft_op = ft_op * padding_op
    vectorial_ft_op = odl.ProductSpaceOperator([[padded_ft_op, 0],
                                                [0, padded_ft_op]])
    ft_kernel = vectorial_ft_op(discretized_kernel)
    return ft_kernel

# Create and initialize deformation field
# Define the momenta and set it to zeroes or ones for test
momenta = vspace.zero()

ft_kernel_fitting = fitting_kernel_2dfft_zero_padding(kernel)
ft_kernel_shape = shape_kernel_2dfft_zero_padding(kernel)

displacement_op = DisplacementOperator(vspace, cptssapce.grid,
                                       discr_space, ft_kernel_fitting)
displ = displacement_op(momenta)


linear_deform_op = LinearizedDeformationOperator(template)
deformed_template = linear_deform_op(displ)

proj_deformed_template = xray_trafo_op(deformed_template)

# Composition of the L2 fitting term with a deformation operator
l2_data_fit_func = L2DataMatchingFunctional(xray_trafo_op.range,
                                            noise_proj_data)

data_fitting_term = l2_data_fit_func * xray_trafo_op * linear_deform_op * displacement_op

# Compute the gradient of shape-based regularization term
kernelmatrix = gaussian_kernel_matrix(cptssapce.grid, sigma)
shape_func = ShapeRegularizationFunctional(vspace, kernelmatrix)
# grad_shape_func = shape_func._gradient(momenta)  # old method
# grad_shape_func = shape_func._gradient_2dfft_zero_padding
# (momenta, kernel)  ## fft method

# Shape regularization parameter, nonnegtive
lambda_shape = 0.0001
# Stepsize for iterations
eta = 50.0
# Iterations for updating alphas
for i in range(2000):
    ft_momenta = vectorial_ft_op(momenta)
    grad_shape_func = vectorial_ft_op_inverse(
        shape_func._gradient_2dfft_zero_padding(ft_momenta, ft_kernel_shape))
    grad_data_fitting_term = data_fitting_term.gradient(momenta)
    momenta -= eta * (
        2 * lambda_shape * grad_shape_func + grad_data_fitting_term)

    if (i+1) % 500 == 0:
        print(i)
#        print(grad_shape_func)
#        print(grad_data_fitting_term)
#        print(momenta)
        displ = displacement_op(momenta)
        print(displ)
        deformed_template = linear_deform_op(displ)  # deformed image
        deformed_template.show(title='Deformed template')
