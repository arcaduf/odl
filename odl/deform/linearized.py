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

"""Linearized deformations based on RKHS vector fields."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

from numbers import Number
import numpy as np

from odl.discr.discr_ops import Gradient
from odl.discr.grid import TensorGrid
from odl.discr.lp_discr import DiscreteLp, DiscreteLpVector
from odl.operator.operator import (
    Operator, OperatorComp, OperatorLeftScalarMult, OperatorRightScalarMult,
    OperatorRightVectorMult,
    _dispatch_call_args, _default_call_in_place, _default_call_out_of_place)
from odl.set.space import LinearSpaceVector
from odl.space.ntuples import Fn
from odl.space.pspace import ProductSpace
from odl.util.utility import preload_first_arg


__all__ = ('Functional', 'FunctionalComp', 'TemplateDeformationOperator',
           'L2DataMatchingFunctional')


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

        x --> I(x + v(x))

    where the vector field ``v`` depends on the deformation parameters
    ``alpha_j`` as follows::

        v(x) = sum_j (K(x, x_j) * alpha_j)

    Here, ``K`` is the RKHS kernel matrix and each ``alpha_j`` is an
    element of ``R^n``.
    """

    def __init__(self, par_space, control_points, template, kernel_func=None,
                 **kwargs):
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
        if isinstance(par_space, Fn):
            # Make a product space with one component
            par_space = ProductSpace(par_space, 1)
        elif isinstance(par_space, ProductSpace):
            pass
        else:
            raise TypeError('expected Rn or ProductSpace as par_space, got '
                            '{!r}.'.format(par_space))

        if not isinstance(template, DiscreteLpVector):
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

        if not isinstance(control_points, TensorGrid):
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
        return isinstance(self.control_points, TensorGrid)

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
            if isinstance(contr_pt_arrs, TensorGrid):
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

    def displacement(self, alphas):
        """Calculate the translation at point x."""
        kernel_mat = self.kernel_matrix()
        return [kernel_mat.T.dot(
                    np.asarray(a).reshape(-1, order=alphas.space[0].order))
                for a in alphas]

    def deform(self, template, alphas):
        """Compute the deformation of the template.

        deformed template: I(x + v(x))

        Parameters
        ----------
        template : `DiscreteLpVector`
            The discrete funtion to be deformed
        alphas: `ProductSpaceVector`
            Deformation parameters for control points. It has ``n``
            components, each of which has size ``N``. Here, ``n`` is
            the number of dimensions and ``N`` the number of control
            points.
        """
        image_pts = self.image_grid.points()
        displacement = self.displacement(alphas)

        for i, v in enumerate(displacement):
            image_pts.T[i] += v

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
        out[:] = self.deform(self.template, alphas)

    def derivative(self, alphas, **kwargs):
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
        deriv_op = TemplateDeformationDerivative(
            alphas, self.control_points, self.template, self.kernel_func,
            cache_kernel_matrix=self._cache_kernel_matrix)

        # Share the cached stuff
        deriv_op._kernel_matrix = self._kernel_matrix
        return deriv_op


class TemplateDeformationDerivative(TemplateDeformationOperator):

    """Frechet derivative of the template deformation operator."""

    def __init__(self, alphas, control_points, template,
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
        """
        super().__init__(alphas.space, control_points, template, kernel_func,
                         **kwargs)
        Operator.__init__(self, alphas.space, template.space, linear=True)
        self.alphas = alphas

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
        grad_space = ProductSpace(self.range, self.ndim)
        if out is None:
            out = grad_space.element()

        for gf, oi in zip(grad_f, out):
            # TODO: do deformation in place
            oi[:] = self.deform(gf, alphas)

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
        grad = Gradient(self.range)
        template_grad = grad(self.template)
        def_grad = self._deform_grad(template_grad, self.alphas,
                                     out=template_grad)

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
            self.alphas, self.control_points, self.template,
            self.kernel_func, cache_kernel_matrix=self._cache_kernel_matrix)

        # Share the cached stuff
        adj_op._kernel_matrix = self._kernel_matrix
        return adj_op


class TemplateDeformationDerivativeAdjoint(TemplateDeformationDerivative):

    """Adjoint of the template deformation operator derivative.

    TODO: write a bit more
    """
    def __init__(self, alphas, control_points, template,
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
        """
        super().__init__(alphas, control_points, template, kernel_func,
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
        grad = Gradient(self.domain)
        template_grad = grad(self.template)
        def_grad = self._deform_grad(template_grad, self.alphas,
                                     out=template_grad)

        for a, gf in zip(out, def_grad):
            self.kernel_matrix().dot(
                b=np.asarray(gf * func).reshape(-1, order=self.domain.order),
                out=np.asarray(a).reshape(-1, order=self.range[0].order))
        out *= self.domain.cell_volume


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
        if not (isinstance(space, DiscreteLp) and space.exponent == 2.0):
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
