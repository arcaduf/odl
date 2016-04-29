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

"""Operators defined on `DiscreteLp`."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

import numpy as np

from odl.discr.lp_discr import DiscreteLp, uniform_discr
from odl.operator.operator import Operator
from odl.space.pspace import ProductSpace


__all__ = ('PartialDerivative', 'Gradient', 'Divergence', 'Laplacian',
           'Resampling', 'ZeroPaddingOperator')


# TODO: make helper function to set edge slices

def finite_diff(f, axis=0, dx=1.0, method='forward', out=None, **kwargs):
    """Calculate the partial derivative of ``f`` along a given ``axis``.

    In the interior of the domain of f, the partial derivative is computed
    using first-order accurate forward or backward difference or
    second-order accurate central differences.

    With padding the same method and thus accuracy is used on endpoints as
    in the interior i.e. forward and backward differences use first-order
    accuracy on edges while central differences use second-order accuracy at
    edges.

    Without padding one-sided forward or backward differences are used at
    the boundaries. The accuracy at the endpoints can then also be
    triggered by the edge order.

    The returned array has the same shape as the input array ``f``.

    Per default forward difference with dx=1 and no padding is used.

    Parameters
    ----------
    f : `array-like`
         An N-dimensional array
    axis : `int`, optional
        The axis along which the partial derivative is evaluated
    dx : `float`, optional
        Scalar specifying the distance between sampling points along ``axis``
    method : {'central', 'forward', 'backward'}, optional
        Finite difference method which is used in the interior of the domain
         of ``f``.
    padding_method : {'constant', 'symmetric'}, optional

        'constant' : Pads values outside the domain of ``f`` with a constant
        value given by ``padding_value``

        'symmetric' : Pads with the reflection of the vector mirrored
        along the edge of the array

        If `None` is given, one-sided forward or backward differences
        are used at the boundary.

    padding_value : `float`, optional
        If ``padding_method`` is 'constant' ``f`` assumes ``padding_value``
        for indices outside the domain of ``f``
    edge_order : {1, 2}, optional
        Edge-order accuracy at the boundaries if no padding is used. If
        `None` the edge-order accuracy at endpoints corresponds to the
        accuracy in the interior. Default: `None`
    out : `numpy.ndarray`, optional
         An N-dimensional array to which the output is written. Has to have
         the same shape as the input array ``f``. Default: `None`

    Returns
    -------
    out : `numpy.ndarray`
        N-dimensional array of the same shape as ``f``. If ``out`` is
        provided, the returned object is a reference to it.

    Notes
    -----
    Without padding the use of second-order accurate edges requires at
    least three elements.

    Central differences with padding cannot be used with first-order
    accurate edges.

    Forward and backward differences with padding use the first-order
    accuracy on edges (as in the interior).

    An edge-order accuracy different from the interior can only be triggered
    without padding i.e. when one-sided differences are used at the edges.

    Examples
    --------
    >>> f = np.array([ 0., 1., 2., 3., 4., 5., 6., 7., 8., 9.])

    >>> finite_diff(f)
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

    Without arguments the above defaults to:

    >>> finite_diff(f, axis=0, dx=1.0, method='forward', padding_method=None,
    ... edge_order=None)
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

    >>> finite_diff(f, dx=0.5)
    array([ 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.])
    >>> finite_diff(f, padding_method='constant')
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., -9.])

    Central differences and different edge orders:

    >>> finite_diff(1/2*f**2, method='central')
    array([-0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.])
    >>> finite_diff(1/2*f**2, method='central', edge_order=1)
    array([ 0.5,  1. ,  2. ,  3. ,  4. ,  5. ,  6. ,  7. ,  8. ,  8.5])

    In-place evaluation:

    >>> out = f.copy()
    >>> out is finite_diff(f, out=out)
    True
    """
    # TODO: implement alternative boundary conditions
    f_arr = np.asarray(f)
    ndim = f_arr.ndim

    if f_arr.shape[axis] < 2:
        raise ValueError('In axis {}: at least two elements required, got {}.'
                         ''.format(axis, f_arr.shape[axis]))

    if axis < 0:
        axis += ndim
    if axis >= ndim:
        raise IndexError('axis {} outside the valid range 0 ... {}'
                         ''.format(axis, ndim - 1))

    if dx <= 0:
        raise ValueError("step length {} not positive.".format(dx))
    else:
        dx = float(dx)

    method, method_in = str(method).lower(), method
    if method not in ('central', 'forward', 'backward'):
        raise ValueError('method {} is not understood'.format(method_in))

    padding_method = kwargs.pop('padding_method', None)
    if padding_method not in ('constant', 'symmetric', None):
        raise ValueError('padding value {} not valid'.format(padding_method))
    if padding_method == 'constant':
        padding_value = float(kwargs.pop('padding_value', 0))

    edge_order = kwargs.pop('edge_order', None)
    if edge_order is None:
        if method == 'central':
            edge_order = 2
        else:
            edge_order = 1
    else:
        if edge_order not in (1, 2):
            raise ValueError('edge order {} not valid'.format(edge_order))

    if out is None:
        out = np.empty_like(f_arr)
    else:
        if out.shape != f.shape:
            raise ValueError('expected output shape {}, got {}.'
                             ''.format(f.shape, out.shape))

    if f_arr.shape[axis] == 2 and edge_order == 2:
        raise ValueError('shape of array to small to use edge order 2')

    if padding_method is not None:
        if method == 'central' and edge_order == 1:
            raise ValueError(
                'central differences with padding cannot be used with '
                'first-order accurate edges')
        if method in ('forward', 'backward') and edge_order == 2:
            raise ValueError('{} differences with padding only use edge '
                             'order 1'.format(method))

    # create slice objects: initially all are [:, :, ..., :]

    # current slice
    slice_out = [slice(None)] * ndim

    # slices used to calculate finite differences
    slice_node1 = [slice(None)] * ndim
    slice_node2 = [slice(None)] * ndim
    slice_node3 = [slice(None)] * ndim

    # Interior of the domain of f

    if method == 'central':
        # 2nd order differences in the interior of the domain of f
        slice_out[axis] = slice(1, -1)
        slice_node1[axis] = slice(2, None)
        slice_node2[axis] = slice(None, -2)
        # 1D equivalent: out[1:-1] = (f[2:] - f[:-2])/2.0
        np.subtract(f_arr[slice_node1], f_arr[slice_node2], out[slice_out])
        out[slice_out] /= 2.0

    elif method == 'forward':
        # 1st order differences in the interior of the domain of f
        slice_out[axis] = slice(1, -1)
        slice_node1[axis] = slice(2, None)
        slice_node2[axis] = slice(1, -1)
        # 1D equivalent: out[1:-1] = (f[2:] - f[1:-1])
        np.subtract(f_arr[slice_node1], f_arr[slice_node2], out[slice_out])

    elif method == 'backward':
        # 1st order differences in the interior of the domain of f
        slice_out[axis] = slice(1, -1)
        slice_node1[axis] = slice(1, -1)
        slice_node2[axis] = slice(None, -2)
        # 1D equivalent: out[1:-1] = (f[1:-1] - f[:-2])
        np.subtract(f_arr[slice_node1], f_arr[slice_node2], out[slice_out])

    # Boundaries

    if padding_method == 'constant':
        # Assume constant value c for indices outside the domain of ``f``

        # With padding the method used on endpoints is the same as in the
        # interior of the domain of f

        if method == 'central':
            # 2nd-order lower edge
            slice_out[axis] = 0
            slice_node1[axis] = 1
            # 1D equivalent: out[0] = (f[1] - c)/2.0
            out[slice_out] = (f_arr[slice_node1] - padding_value) / 2.0

            # 2nd-order upper edge
            slice_out[axis] = -1
            slice_node2[axis] = -2
            # 1D equivalent: out[-1] = (c - f[-2])/2.0
            out[slice_out] = (padding_value - f_arr[slice_node2]) / 2.0

        elif method == 'forward':
            # 1st-oder lower edge
            slice_out[axis] = 0
            slice_node1[axis] = 1
            slice_node2[axis] = 0
            # 1D equivalent: out[0] = f[1] - f[0]
            out[slice_out] = f_arr[slice_node1] - f_arr[slice_node2]

            # 1st-oder upper edge
            slice_out[axis] = -1
            slice_node2[axis] = -1
            # 1D equivalent: out[-1] = c - f[-1]
            out[slice_out] = padding_value - f_arr[slice_node2]

        elif method == 'backward':
            # 1st-oder lower edge
            slice_out[axis] = 0
            slice_node1[axis] = 0
            # 1D equivalent: out[0] = f[0] - c
            out[slice_out] = f_arr[slice_node1] - padding_value

            # 1st-oder upper edge
            slice_out[axis] = -1
            slice_node1[axis] = -1
            slice_node2[axis] = -2
            # 1D equivalent: out[-1] = f[-1] - f[-2]
            out[slice_out] = f_arr[slice_node1] - f_arr[slice_node2]

    elif padding_method == 'symmetric':
        # Values of f for indices outside the domain of f are replicates of
        # the edge values

        # With padding the method used on endpoints is the same as in the
        # interior of the domain of f

        if method == 'central':
            # 2nd-order lower edge
            slice_out[axis] = 0
            slice_node1[axis] = 1
            slice_node2[axis] = 0
            # 1D equivalent: out[0] = (f[1] - f[0])/2.0
            out[slice_out] = (f_arr[slice_node1] - f_arr[slice_node2]) / 2.0

            # 2nd-order upper edge
            slice_out[axis] = -1
            slice_node1[axis] = -1
            slice_node2[axis] = -2
            # 1D equivalent: out[-1] = (f[-1] - f[-2])/2.0
            out[slice_out] = (f_arr[slice_node1] - f_arr[slice_node2]) / 2.0

        elif method == 'forward':
            # 1st-oder lower edge
            slice_out[axis] = 0
            slice_node1[axis] = 1
            slice_node2[axis] = 0
            # 1D equivalent: out[0] = f[1] - f[0]
            out[slice_out] = f_arr[slice_node1] - f_arr[slice_node2]

            # 1st-oder upper edge
            slice_out[axis] = -1
            # 1D equivalent: out[-1] = f[-1] - f[-1] = 0
            out[slice_out] = 0

        elif method == 'backward':
            # 1st-oder lower edge
            slice_out[axis] = 0
            # 1D equivalent: out[0] = f[0] - f[0] = 0
            out[slice_out] = 0

            # 1st-oder upper edge
            slice_out[axis] = -1
            slice_node1[axis] = -1
            slice_node2[axis] = -2
            # 1D equivalent: out[-1] = f[-1] - f[-2]
            out[slice_out] = f_arr[slice_node1] - f_arr[slice_node2]

    # Use one-sided differences on the endpoints
    else:

        # Edge-order accuracy is triggered implicitly by the method used or
        # explicitly using ``edge_order``

        # 1st order edges
        if edge_order == 1:
            # lower boundary
            slice_out[axis] = 0
            slice_node1[axis] = 1
            slice_node2[axis] = 0
            # 1D equivalent: out[0] = (f[1] - f[0])
            out[slice_out] = f_arr[slice_node1] - f_arr[slice_node2]

            # upper boundary
            slice_out[axis] = -1
            slice_node1[axis] = -1
            slice_node2[axis] = -2
            # 1D equivalent: out[-1] = (f[-1] - f[-2])
            out[slice_out] = f_arr[slice_node1] - f_arr[slice_node2]

        # 2nd order edges
        elif edge_order == 2:
            # lower boundary
            slice_out[axis] = 0
            slice_node1[axis] = 0
            slice_node2[axis] = 1
            slice_node3[axis] = 2
            # 1D equivalent: out[0] = -(3*f[0] - 4*f[1] + f[2]) / 2.0
            out[slice_out] = -(3.0 * f_arr[slice_node1] - 4.0 * f_arr[
                slice_node2] + f_arr[slice_node3]) / 2.0

            # upper boundary
            slice_out[axis] = -1
            slice_node1[axis] = -1
            slice_node2[axis] = -2
            slice_node3[axis] = -3
            # 1D equivalent: out[-1] = (3*f[-1] - 4*f[-2] + f[-3]) / 2.0
            out[slice_out] = (3.0 * f_arr[slice_node1] - 4.0 * f_arr[
                slice_node2] + f_arr[slice_node3]) / 2.0

    # divide by step size
    out /= dx

    return out


class PartialDerivative(Operator):

    """Calculate the discrete partial derivative along a given axis.

    Calls helper function `finite_diff` to calculate finite difference.
    Preserves the shape of the underlying grid.
    """
    # TODO: implement adjoint

    def __init__(self, space, axis=0, method='forward', padding_method=None,
                 padding_value=0, edge_order=None):
        """Initialize an operator instance.

        Parameters
        ----------
        space : `DiscreteLp`
            The space of elements which the operator is acting on
        axis : `int`, optional
            The axis along which the partial derivative is evaluated
        method : {'central', 'forward', 'backward'}, optional
            Finite difference method which is used in the interior of the
            domain of ``f``
        padding_method : {'constant', 'symmetric'}, optional

            'constant' : Pads values outside the domain of ``f`` with a
            constant value given by ``padding_value``

            'symmetric' : Pads with the reflection of the vector mirrored
            along the edge of the array

            If `None` is given, one-sided forward or backward differences
            are used at the boundary

        padding_value : `float`, optional
            If ``padding_method`` is 'constant' ``f`` assumes
            ``padding_value`` for indices outside the domain of ``f``
        edge_order : {1, 2}, optional
            Edge-order accuracy at the boundaries if no padding is used. If
            `None` the edge-order accuracy at endpoints corresponds to the
            accuracy in the interior.
        """
        if not isinstance(space, DiscreteLp):
            raise TypeError('space {!r} is not a DiscreteLp '
                            'instance.'.format(space))

        super().__init__(domain=space, range=space, linear=True)
        self.axis = axis
        self.dx = space.cell_sides[axis]
        self.method = method
        self.padding_method = padding_method
        self.padding_value = padding_value
        self.edge_order = edge_order

    def _call(self, x, out=None):
        """Apply gradient operator to ``x`` and store result in ``out``.

        Parameters
        ----------
        x : ``domain`` `element`
            Input vector to which the operator is applied to
        out : ``range`` element, optional
            Output vector to which the result is written

        Returns
        -------
        out : ``range`` `element`
            Result of the evaluation. If ``out`` is provided, the
            returned object is a reference to it.

        Examples
        --------
        >>> from odl import uniform_discr
        >>> data = np.array([[ 0.,  1.,  2.,  3.,  4.],
        ...                  [ 0.,  2.,  4.,  6.,  8.]])
        >>> discr = uniform_discr([0, 0], [2, 1], data.shape)
        >>> par_deriv = PartialDerivative(discr)
        >>> f = par_deriv.domain.element(data)
        >>> par_div_f = par_deriv(f)
        >>> print(par_div_f)
        [[0.0, 1.0, 2.0, 3.0, 4.0],
         [0.0, 1.0, 2.0, 3.0, 4.0]]
        """
        if out is None:
            out = self.range.element()

        # TODO: this pipes CUDA arrays through NumPy. Write native operator.
        out_arr = out.asarray()
        finite_diff(x.asarray(), out=out_arr, axis=self.axis, dx=self.dx,
                    method=self.method, padding_method=self.padding_method,
                    padding_value=self.padding_value,
                    edge_order=self.edge_order)

        # self assignment: no overhead in the case asarray is a view
        out[:] = out_arr
        return out

    @property
    def adjoint(self):
        """Return the adjoint operator."""
        raise NotImplementedError('adjoint not implemented')


class Gradient(Operator):
    """Spatial gradient operator for `DiscreteLp` spaces.

    Calls helper function `finite_diff` to calculate each component of the
    resulting product space vector. For the adjoint of the `Gradient`
    operator, zero padding is assumed to match the negative `Divergence`
    operator
    """

    def __init__(self, space, method='forward'):
        """Initialize a `Gradient` operator instance.

        Zero padding is assumed for the adjoint of the `Gradient`
        operator to match  negative `Divergence` operator.

        Parameters
        ----------
        space : `DiscreteLp`
            The space of elements which the operator is acting on.
        method : {'central', 'forward', 'backward'}, optional
            Finite difference method to be used
        """

        if not isinstance(space, DiscreteLp):
            raise TypeError('space {!r} is not a `DiscreteLp` '
                            'instance.'.format(space))

        self.method = method
        super().__init__(
            domain=space, range=ProductSpace(space, space.ndim), linear=True)

    def _call(self, x, out=None):
        """Calculate the spatial gradient of ``x``.

        Parameters
        ----------
        x : ``domain`` `element`
            Input vector to which the `Gradient` operator is applied
        out : ``range`` `element`, optional
            Output vector to which the result is written

        Returns
        -------
        out : ``range`` `element`
            Result of the evaluation. If ``out`` is provided, the returned
            object is a reference to it.

        Examples
        --------
        >>> from odl import uniform_discr
        >>> data = np.array([[ 0., 1., 2., 3., 4.],
        ...                  [ 0., 2., 4., 6., 8.]])
        >>> discr = uniform_discr([0, 0], [2, 5], data.shape)
        >>> f = discr.element(data)
        >>> grad = Gradient(discr)
        >>> grad_f = grad(f)
        >>> print(grad_f[0])
        [[0.0, 1.0, 2.0, 3.0, 4.0],
         [0.0, -2.0, -4.0, -6.0, -8.0]]
        >>> print(grad_f[1])
        [[1.0, 1.0, 1.0, 1.0, -4.0],
         [2.0, 2.0, 2.0, 2.0, -8.0]]

        Verify adjoint:

        >>> g = grad.range.element((data, data ** 2))
        >>> adj_g = grad.adjoint(g)
        >>> print(adj_g)
        [[0.0, -2.0, -5.0, -8.0, -11.0],
         [0.0, -5.0, -14.0, -23.0, -32.0]]
        >>> g.inner(grad_f) / f.inner(adj_g)
        1.0
        """
        if out is None:
            out = self.range.element()

        x_arr = x.asarray()
        ndim = self.domain.ndim
        dx = self.domain.cell_sides

        for axis in range(ndim):
            out_arr = out[axis].asarray()

            finite_diff(x_arr, axis=axis, dx=dx[axis], method=self.method,
                        padding_method='constant', padding_value=0,
                        out=out_arr, )

            out[axis][:] = out_arr

        return out

    @property
    def adjoint(self):
        """Return the adjoint operator.

        Assuming implicit zero padding, the adjoint operator is given by the
        negative of the `Divergence` operator.

        The Divergence is constructed from a ``space`` as a product space
        operator ``space^n --> space``, hence we need to provide the domain of
        this operator.
        """
        if self.method == 'central':
            return - Divergence(self.domain, 'central')
        elif self.method == 'forward':
            return - Divergence(self.domain, 'backward')
        elif self.method == 'backward':
            return - Divergence(self.domain, 'forward')
        else:
            return super().adjoint


class Divergence(Operator):
    """Divergence operator for `DiscreteLp` spaces.

    Calls helper function `finite_diff` for each component of the input
    product space vector. For the adjoint of the `Divergence` operator to
    match the negative `Gradient` operator implicit zero is assumed.
    """

    def __init__(self, space, method='forward'):
        """Initialize a `Divergence` operator instance.

        Zero padding is assumed for the adjoint of the `Divergence`
        operator to match negative `Gradient` operator.

        Parameters
        ----------
        space : `DiscreteLp`
            The space of elements which the operator is acting on
        method : {'central', 'forward', 'backward'}, optional
            Finite difference method to be used
        """
        if not isinstance(space, DiscreteLp):
            raise TypeError('space {!r} is not a `DiscreteLp` '
                            'instance.'.format(space))

        self.method = method
        super().__init__(domain=ProductSpace(space, space.ndim),
                         range=space, linear=True)

    def _call(self, x, out=None):
        """Calculate the divergence of ``x``.

        Parameters
        ----------
        x : ``domain`` `element`
            `ProductSpaceVector` to which the divergence operator
            is applied
        out : ``range`` `element`, optional
            Output vector to which the result is written

        Returns
        -------
        out : ``range`` `element`
            Result of the evaluation. If ``out`` is provided, the returned
            object is a reference to it.

        Examples
        --------
        >>> from odl import uniform_discr
        >>> data = np.array([[0., 1., 2., 3., 4.],
        ...                  [1., 2., 3., 4., 5.],
        ...                  [2., 3., 4., 5., 6.]])
        >>> space = uniform_discr([0, 0], [3, 5], data.shape)
        >>> div = Divergence(space)
        >>> f = div.domain.element([data, data])
        >>> div_f = div(f)
        >>> print(div_f)
        [[2.0, 2.0, 2.0, 2.0, -3.0],
         [2.0, 2.0, 2.0, 2.0, -4.0],
         [-1.0, -2.0, -3.0, -4.0, -12.0]]

        Verify adjoint:

        >>> g = div.range.element(data ** 2)
        >>> adj_div_g = div.adjoint(g)
        >>> g.inner(div_f) / f.inner(adj_div_g)
        1.0
        """
        if out is None:
            out = self.range.element()

        ndim = self.range.ndim
        dx = self.range.cell_sides

        out_arr = out.asarray()
        tmp = np.empty(out.shape, out.dtype, order=out.space.order)
        for axis in range(ndim):
            finite_diff(x[axis], axis=axis, dx=dx[axis], method=self.method,
                        padding_method='constant', padding_value=0, out=tmp)
            if axis == 0:
                out_arr[:] = tmp
            else:
                out_arr += tmp

        # self assignment: no overhead in the case asarray is a view
        out[:] = out_arr
        return out

    @property
    def adjoint(self):
        """Return the adjoint operator.

        Assuming implicit zero padding, the adjoint operator is given by the
        negative of the `Gradient` operator.
        """
        if self.method == 'central':
            return - Gradient(self.range, 'central')
        elif self.method == 'forward':
            return - Gradient(self.range, 'backward')
        elif self.method == 'backward':
            return - Gradient(self.range, 'forward')
        else:
            return super().adjoint


class Laplacian(Operator):
    """Spatial Laplacian operator for `DiscreteLp` spaces.

    Calls helper function `finite_diff` to calculate each component of the
    resulting product space vector.

    Outside the domain zero padding is assumed.
    """

    def __init__(self, space):
        """Initialize a `Laplacian` operator instance.

        Parameters
        ----------
        space : `DiscreteLp`
            The space of elements which the operator is acting on
        """

        if not isinstance(space, DiscreteLp):
            raise TypeError('space {!r} is not a `DiscreteLp` '
                            'instance.'.format(space))

        super().__init__(domain=space, range=space, linear=True)

    def _call(self, x, out=None):
        """Calculate the spatial Laplacian of ``x``.

        Parameters
        ----------
        x : ``domain`` `element`
            Input vector to which the `Laplacian` operator is
            applied
        out : ``range`` `element`, optional
            Output vector to which the result is written

        Returns
        -------
        out : ``range`` `element`
            Result of the evaluation. If ``out`` is provided, the returned
            object is a reference to it.

        Examples
        --------
        >>> from odl import uniform_discr
        >>> data = np.array([[ 0., 0., 0.],
        ...                  [ 0., 1., 0.],
        ...                  [ 0., 0., 0.]])
        >>> space = uniform_discr([0, 0], [3, 3], data.shape)
        >>> f = space.element(data)
        >>> lap = Laplacian(space)
        >>> print(lap(f))
        [[0.0, 1.0, 0.0],
         [1.0, -4.0, 1.0],
         [0.0, 1.0, 0.0]]
        """
        if out is None:
            out = self.range.zero()
        else:
            out.set_zero()

        x_arr = x.asarray()
        out_arr = out.asarray()
        tmp = np.empty(out.shape, out.dtype, order=out.space.order)

        ndim = self.domain.ndim
        dx = self.domain.cell_sides

        for axis in range(ndim):
            # TODO: this can be optimized

            finite_diff(x_arr, axis=axis, dx=dx[axis] ** 2,
                        method='forward', padding_method='constant',
                        padding_value=0, out=tmp)

            out_arr[:] += tmp

            finite_diff(x_arr, axis=axis, dx=dx[axis] ** 2,
                        method='backward', padding_method='constant',
                        padding_value=0, out=tmp)

            out_arr[:] -= tmp

        out[:] = out_arr
        return out

    @property
    def adjoint(self):
        """Return the adjoint operator.

        The laplacian is self-adjoint, so this returns ``self``.
        """
        return self


class Resampling(Operator):
    """A operator that resamples a vector on another grid.

    The operator uses the underlying `Discretization.sampling` and
    `Discretization.interpolation` operators to achieve this.

    The spaces need to have the same `Discretization.uspace` in order for this
    to work. The dspace types may however be different, although performance
    may vary drastically.
    """

    def __init__(self, domain, range):
        """Initialize a Resampling.

        Parameters
        ----------
        domain : `LinearSpace`
            The space that should be cast from
        range : `LinearSpace`
            The space that should be cast to

        Examples
        --------
        Create two spaces with different number of points and create resampling
        operator.

        >>> import odl
        >>> X = odl.uniform_discr(0, 1, 3)
        >>> Y = odl.uniform_discr(0, 1, 6)
        >>> resampling = Resampling(X, Y)
        """
        if domain.uspace != range.uspace:
            raise ValueError('domain.uspace ({}) does not match range.uspace '
                             '({})'.format(domain.uspace, range.uspace))

        super().__init__(domain=domain, range=range, linear=True)

    def _call(self, x, out=None):
        """Apply resampling operator.

        The vector ``x`` is resampled using the sampling and interpolation
        operators of the underlying spaces.

        Examples
        --------
        Create two spaces with different number of points and create resampling
        operator. Apply operator to vector.

        >>> import odl
        >>> X = odl.uniform_discr(0, 1, 3)
        >>> Y = odl.uniform_discr(0, 1, 6)
        >>> resampling = Resampling(X, Y)
        >>> print(resampling([0, 1, 0]))
        [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]

        The result depends on the interpolation chosen for the underlying
        spaces.

        >>> Z = odl.uniform_discr(0, 1, 3, interp='linear')
        >>> linear_resampling = Resampling(Z, Y)
        >>> print(linear_resampling([0, 1, 0]))
        [0.0, 0.25, 0.75, 0.75, 0.25, 0.0]
        """
        if out is None:
            return x.interpolation
        else:
            out.sampling(x.interpolation)

    @property
    def inverse(self):
        """Return an (approximate) inverse.

        Returns
        -------
        inverse : Resampling
            The resampling operator defined in the inverse direction.

        See Also
        --------
        adjoint : resampling is unitary, so adjoint is inverse.
        """
        return Resampling(self.range, self.domain)

    @property
    def adjoint(self):
        """Return an (approximate) adjoint.

        The result is only exact if the interpolation and sampling operators
        of the underlying spaces match exactly.

        Returns
        -------
        adjoint : Resampling
            The resampling operator defined in the inverse direction.

        Examples
        --------
        Create resampling operator and inverse

        >>> import odl
        >>> X = odl.uniform_discr(0, 1, 3)
        >>> Y = odl.uniform_discr(0, 1, 6)
        >>> resampling = Resampling(X, Y)
        >>> resampling_inv = resampling.inverse

        The inverse is proper left inverse if the resampling goes from a
        lower sampling to a higher sampling

        >>> x = [0.0, 1.0, 0.0]
        >>> print(resampling_inv(resampling(x)))
        [0.0, 1.0, 0.0]

        But can fail in the other direction

        >>> y = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        >>> print(resampling(resampling_inv(y)))
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        """
        return self.inverse

class ZeroPaddingOperator(Operator):

    def __init__(self, space, padfrac, axes=None):
        assert isinstance(space, DiscreteLp)
        assert space.partition.is_regular
        assert np.all(np.greater_equal(padfrac, 0))

        if axes is None:
            axes = np.arange(space.ndim)

        self._padfrac = padfrac
        self._axes = axes
        padding = []
        for i, p in enumerate(padfrac):
            if i in axes:
                padding.append(p)
            else:
                padding.append(0.0)

        padding = np.array(padding, dtype=float)

        n_old = np.array(space.shape, dtype=int)
        n_new = np.ceil((1 + padding) * n_old).astype(int)
        n_diff = n_new - n_old
        n_left = n_diff // 2

        self._slc_start = n_left
        self._slc_end = self._slc_start + n_old

        stride = space.cell_sides
        newbegin = space.domain.begin - n_left * stride
        newend = space.domain.end + (n_diff - n_left) * stride

        # TODO: take care of boundary nodes
        # TODO: handle CUDA
        extended_space = uniform_discr(
            min_corner=newbegin, max_corner=newend, nsamples=n_new,
            exponent=space.exponent, interp=space.interp, impl='numpy')

        super().__init__(domain=space, range=extended_space, linear=True)

    def _call(self, f):
        nstart = self._slc_start
        nend = self._slc_end
        slc = []
        for s, e in zip(nstart, nend):
            slc.append(slice(s, e))

        padded_arr = self.range.zero().asarray()
        padded_arr[slc] = f
        return padded_arr

    @property
    def inverse(self):

        zeropad_op = self

        class _Inverse(ZeroPaddingOperator):
            def __init__(self):
                super().__init__(zeropad_op.domain, zeropad_op._padfrac,
                                 axes=zeropad_op._axes)
                self._domain, self._range = self._range, self._domain

            def _call(self, f):
                nstart = self._slc_start
                nend = self._slc_end
                slc = []
                for s, e in zip(nstart, nend):
                    slc.append(slice(s, e))

                return f.asarray()[slc]

        return _Inverse()

    @property
    def adjoint(self):
        if self.domain.exponent != 2.0:
            raise NotImplementedError
        return self.inverse

if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
