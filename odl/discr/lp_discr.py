﻿# Copyright 2014-2016 The ODL development group
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

""":math:`L^p` type discretizations of function spaces."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super, str

import numpy as np
from numbers import Integral

from odl.discr.discretization import (
    DiscretizedSpace, DiscretizedSpaceVector, dspace_type)
from odl.discr.discr_mappings import (
    PointCollocation, NearestInterpolation, LinearInterpolation,
    PerAxisInterpolation)
from odl.discr.partition import (
    RectPartition, uniform_partition_fromintv, uniform_partition)
from odl.set.sets import RealNumbers, ComplexNumbers
from odl.set.domain import IntervalProd
from odl.space.cu_ntuples import CUDA_AVAILABLE, CudaFn
from odl.space.fspace import FunctionSpace
from odl.space.ntuples import Fn
from odl.util.normalize import (
    normalized_scalar_param_list, safe_int_conv, normalized_nodes_on_bdry)
from odl.util.numerics import apply_on_boundary
from odl.util.ufuncs import DiscreteLpUFuncs
from odl.util.utility import (
    is_real_dtype, is_complex_floating_dtype, dtype_repr)

__all__ = ('DiscreteLp', 'DiscreteLpVector',
           'uniform_discr_frompartition', 'uniform_discr_fromspace',
           'uniform_discr_fromintv', 'uniform_discr',
           'uniform_discr_fromdiscr', 'discr_sequence_space')

_SUPPORTED_INTERP = ('nearest', 'linear')


class DiscreteLp(DiscretizedSpace):

    """Discretization of a Lebesgue :math:`L^p` space."""

    def __init__(self, fspace, partition, dspace, exponent=2.0,
                 interp='nearest', **kwargs):
        """Initialize a new instance.

        Parameters
        ----------
        fspace : `FunctionSpace`
            The continuous space to be discretized
        partition : `RectPartition`
            Partition of (a subset of) ``fspace.domain`` based on a
            `TensorGrid`
        dspace : `FnBase`
            Space of elements used for data storage. It must have the
            same `FnBase.field` as ``fspace``
        exponent : positive `float`, optional
            The parameter :math:`p` in :math:`L^p`. If the exponent is
            not equal to the default 2.0, the space has no inner
            product.
        interp : `str` or `sequence` of `str`, optional
            The interpolation type to be used for discretization.
            A sequence is interpreted as interpolation scheme per
            axis.

            'nearest' : use nearest-neighbor interpolation (default)

            'linear' : use linear interpolation

        order : {'C', 'F'}, optional
            Ordering of the axes in the data storage. 'C' means the
            first axis varies slowest, the last axis fastest;
            vice versa for 'F'.
            Default: 'C'
        """
        if not isinstance(fspace, FunctionSpace):
            raise TypeError('{!r} is not a FunctionSpace instance'
                            ''.format(fspace))
        if not isinstance(fspace.domain, IntervalProd):
            raise TypeError('function space domain {!r} is not an '
                            'IntervalProd instance'.format(fspace.domain))
        if not isinstance(partition, RectPartition):
            raise TypeError('`partition` {!r} is not a RectPartition '
                            'instance'.format(partition))
        if not fspace.domain.contains_set(partition.set):
            raise ValueError('`partition` {} is not a subset of the function '
                             'domain {}'.format(partition, fspace.domain))
        if fspace.out_dtype != dspace.dtype:
            raise ValueError('`out_dtype` of Function Space {} does not match '
                             'dspace dtype {}'
                             ''.format(fspace.out_dtype, dspace.dtype))

        try:
            # Got single string
            interp, interp_in = str(interp + '').lower(), interp
            if interp not in _SUPPORTED_INTERP:
                raise ValueError("`interp` type '{}' not understood"
                                 "".format(interp_in))
            self._interp = [interp] * partition.ndim
        except TypeError:
            # Got sequence of strings
            if len(interp) != partition.ndim:
                raise ValueError('expected {} (ndim) entries in interp, '
                                 'got {}'.format(partition.ndim, len(interp)))

            self._interp = [str(s).lower() for s in interp]
            if any(s not in _SUPPORTED_INTERP for s in self._interp):
                raise ValueError('interp sequence {} contains illegal '
                                 'values'.format(interp))

        order = str(kwargs.pop('order', 'C'))
        if str(order).upper() not in ('C', 'F'):
            raise ValueError('`order` {!r} not recognized'.format(order))
        else:
            self._order = str(order).upper()

        self._partition = partition
        sampling = PointCollocation(fspace, self.partition, dspace,
                                    order=self.order)
        if all(s == 'nearest' for s in self.interp):
            interpol = NearestInterpolation(fspace, self.partition, dspace,
                                            order=self.order)
        elif all(s == 'linear' for s in self.interp):
            interpol = LinearInterpolation(fspace, self.partition, dspace,
                                           order=self.order)
        else:
            interpol = PerAxisInterpolation(
                fspace, self.partition, dspace, self.interp, order=self.order)

        DiscretizedSpace.__init__(self, fspace, dspace, sampling, interpol)
        self._exponent = float(exponent)
        if (hasattr(self.dspace, 'exponent') and
                self.exponent != dspace.exponent):
            raise ValueError('`exponent` {} not equal to data space exponent '
                             '{}'.format(self.exponent, dspace.exponent))

        if self.field == RealNumbers():
            self._real_space = self
            self._complex_space = None
        elif self.field == ComplexNumbers():
            self._real_space = None
            self._complex_space = self
        else:
            self._real_space = None
            self._complex_space = None

    @property
    def min_corner(self):
        """Vector of minimal coordinates of the function domain."""
        return self.partition.begin

    @property
    def max_corner(self):
        """Vector of maximal coordinates of the function domain."""
        return self.partition.end

    @property
    def order(self):
        """Axis ordering for array flattening."""
        return self._order

    @property
    def partition(self):
        """The `RectPartition` of the domain."""
        return self._partition

    @property
    def is_uniform(self):
        """Return `True` if ``self.partition`` is uniform."""
        return self.partition.is_uniform

    @property
    def grid(self):
        """Sampling grid of the discretization mappings."""
        return self.partition.grid

    @property
    def shape(self):
        """Shape of the underlying partition."""
        return self.partition.shape

    @property
    def ndim(self):
        """Number of dimensions."""
        return self.partition.ndim

    @property
    def size(self):
        """Total number of underlying partition cells."""
        return self.partition.size

    @property
    def cell_sides(self):
        """Side lengths of a cell in an underlying *uniform* partition."""
        return self.partition.cell_sides

    @property
    def cell_volume(self):
        """Cell volume of an underlying regular partition."""
        return self.partition.cell_volume

    @property
    def meshgrid(self):
        """All sampling points in the partition as a sparse meshgrid."""
        return self.partition.meshgrid

    def points(self):
        """All sampling points in the partition."""
        return self.partition.points()

    @property
    def exponent(self):
        """The exponent ``p`` in ``L^p``."""
        return self._exponent

    def element(self, inp=None, **kwargs):
        """Create an element from ``inp`` or from scratch.

        Parameters
        ----------
        inp : optional
            The input data to create an element from. It needs to be
            understood by either the `sampling` operator of this
            instance or by its ``dspace.element`` method.
        kwargs :
            Additional arguments passed on to `sampling` when called
            on ``inp``, in the form ``sampling(inp, **kwargs)``.
            This can be used e.g. for functions with parameters.

        Returns
        -------
        element : `DiscretizedSetVector`
            The discretized element, calculated as ``sampling(inp)`` or
            ``dspace.element(inp)``, tried in this order.

        Examples
        --------
        Elements can be created from array-like objects that represent
        an already discretized function:

        >>> import odl
        >>> space = odl.uniform_discr(-1, 1, 4)
        >>> space.element([1, 2, 3, 4])
        uniform_discr(-1.0, 1.0, 4).element([1.0, 2.0, 3.0, 4.0])
        >>> vector = odl.Rn(4).element([0, 1, 2, 3])
        >>> space.element(vector)
        uniform_discr(-1.0, 1.0, 4).element([0.0, 1.0, 2.0, 3.0])

        On the other hand, non-discretized objects like Python functions
        can be discretized "on the fly":

        >>> def f(x):
        ...     return x * 2
        ...
        >>> space.element(f)
        uniform_discr(-1.0, 1.0, 4).element([-1.5, -0.5, 0.5, 1.5])

        This works also with parameterized functions, however only
        through keyword arguments (not positional arguments with
        defaults):

        >>> def f(x, **kwargs):
        ...     c = kwargs.pop('c', 0.0)
        ...     return np.maximum(x, c)
        ...
        >>> space = odl.uniform_discr(-1, 1, 4)
        >>> space.element(f, c=0.5)
        uniform_discr(-1.0, 1.0, 4).element([0.5, 0.5, 0.5, 0.75])

        See also
        --------
        sampling : create a discrete element from an undiscretized one
        """
        if inp is None:
            return self.element_type(self, self.dspace.element())
        elif inp in self:
            return inp
        elif inp in self.dspace:
            return self.element_type(self, inp)
        elif callable(inp):
            # uspace element -> discretize
            inp_elem = self.uspace.element(inp)
            return self.element_type(self, self.sampling(inp_elem, **kwargs))
        else:
            # Sequence-type input
            arr = np.asarray(inp, dtype=self.dtype, order=self.order)
            input_shape = arr.shape

            # Attempt to handle some cases with exact match
            if arr.ndim > 1 and self.ndim == 1:
                arr = np.squeeze(arr)  # Squeeze could solve the problem
                if arr.shape != self.shape:
                    raise ValueError(
                        'input shape {} does not match grid shape {}'
                        ''.format(input_shape, self.shape))
            elif arr.shape != self.shape and arr.ndim == self.ndim:
                # Try to broadcast the array if possible
                try:
                    arr = np.broadcast_to(arr, self.shape)
                except AttributeError:
                    # The above requires numpy 1.10, fallback impl else
                    shape = [m if n == 1 and m != 1 else 1
                             for n, m in zip(arr.shape, self.shape)]
                    arr = arr + np.zeros(shape, dtype=arr.dtype)
            arr = arr.ravel(order=self.order)
            return self.element_type(self, self.dspace.element(arr))

    @property
    def interp(self):
        """Interpolation type of this discretization."""
        return self._interp

    def _astype(self, dtype):
        """Internal helper for ``astype``."""
        fspace = self.uspace.astype(dtype)
        dspace = self.dspace.astype(dtype)
        return type(self)(fspace, self.partition, dspace,
                          exponent=self.exponent, interp=self.interp,
                          order=self.order)

    # Overrides for space functions depending on partition
    #
    # The inherited methods by default use a weighting by a constant
    # (the grid cell size). In dimensions where the partitioned set contains
    # only a fraction of the outermost cells (e.g. if the outermost grid
    # points lie at the boundary), the corresponding contribuitons to
    # discretized integrals need to be scaled by that fraction.
    def _inner(self, x, y):
        """Return ``self.inner(x, y)``."""
        bdry_fracs = self.partition.boundary_cell_fractions
        if (np.allclose(bdry_fracs, 1.0) or
                self.exponent == float('inf') or
                not getattr(self.dspace, 'is_weighted', False)):
            # no boundary weighting
            return super()._inner(x, y)
        else:
            # TODO: implement without copying x
            func_list = _scaling_func_list(bdry_fracs)

            x_arr = apply_on_boundary(x, func=func_list, only_once=False)
            return super()._inner(self.element(x_arr), y)

    def _norm(self, x):
        """Return ``self.norm(x)``."""
        bdry_fracs = self.partition.boundary_cell_fractions
        if (np.allclose(bdry_fracs, 1.0) or
                self.exponent == float('inf') or
                not getattr(self.dspace, 'is_weighted', False)):
            # no boundary weighting
            return super()._norm(x)
        else:
            # TODO: implement without copying x
            func_list = _scaling_func_list(bdry_fracs, exponent=self.exponent)

            x_arr = apply_on_boundary(x, func=func_list, only_once=False)
            return super()._norm(self.element(x_arr))

    def _dist(self, x, y):
        """Return ``self.dist(x, y)``."""
        bdry_fracs = self.partition.boundary_cell_fractions
        if (np.allclose(bdry_fracs, 1.0) or
                self.exponent == float('inf') or
                not getattr(self.dspace, 'is_weighted', False)):
            # no boundary weighting
            return super()._dist(x, y)
        else:
            # TODO: implement without copying x
            func_list = _scaling_func_list(bdry_fracs, exponent=self.exponent)

            arrs = [apply_on_boundary(vec, func=func_list, only_once=False)
                    for vec in (x, y)]

            return super()._dist(self.element(arrs[0]), self.element(arrs[1]))

    def __repr__(self):
        """Return ``repr(self)``."""
        # Check if the factory repr can be used
        if (uniform_partition_fromintv(
                self.uspace.domain, self.shape,
                nodes_on_bdry=False) == self.partition):
            use_uniform = True
            nodes_on_bdry = False
        elif (uniform_partition_fromintv(
                self.uspace.domain, self.shape,
                nodes_on_bdry=True) == self.partition):
            use_uniform = True
            nodes_on_bdry = True
        else:
            use_uniform = False

        if use_uniform:
            if isinstance(self.dspace, Fn):
                impl = 'numpy'
                default_dtype = np.float64
            elif isinstance(self.dspace, CudaFn):
                impl = 'cuda'
                default_dtype = np.float32
            else:  # This should never happen
                raise RuntimeError('unable to determine data space impl')
            arg_fstr = '{}, {}, {}'
            if self.exponent != 2.0:
                arg_fstr += ', exponent={exponent}'
            if self.dtype != default_dtype:
                arg_fstr += ', dtype={dtype}'
            if not all(s == 'nearest' for s in self.interp):
                arg_fstr += ', interp={interp!r}'
            if impl != 'numpy':
                arg_fstr += ', impl={impl!r}'
            if self.order != 'C':
                arg_fstr += ', order={order!r}'
            if nodes_on_bdry:
                arg_fstr += ', nodes_on_bdry={nodes_on_bdry!r}'

            if self.ndim == 1:
                min_str = '{!r}'.format(self.uspace.domain.min()[0])
                max_str = '{!r}'.format(self.uspace.domain.max()[0])
                shape_str = '{!r}'.format(self.shape[0])
            else:
                min_str = '{!r}'.format(list(self.uspace.domain.min()))
                max_str = '{!r}'.format(list(self.uspace.domain.max()))
                shape_str = '{!r}'.format(list(self.shape))

            arg_str = arg_fstr.format(
                min_str, max_str, shape_str,
                exponent=self.exponent,
                dtype=dtype_repr(self.dtype),
                interp=self.interp,
                impl=impl,
                order=self.order,
                nodes_on_bdry=nodes_on_bdry)
            return 'uniform_discr({})'.format(arg_str)
        else:
            arg_fstr = '''
    {!r},
    {!r},
    {!r}
    '''
            if self.exponent != 2.0:
                arg_fstr += ', exponent={ex}'
            if self.interp != 'nearest':
                arg_fstr += ', interp={interp!r}'
            if self.order != 'C':
                arg_fstr += ', order={order!r}'

            arg_str = arg_fstr.format(
                self.uspace, self.partition, self.dspace, interp=self.interp,
                ex=self.exponent, order=self.order)
            return '{}({})'.format(self.__class__.__name__, arg_str)

    def __str__(self):
        """Return ``str(self)``."""
        return self.__repr__()

    @property
    def element_type(self):
        """`DiscreteLpVector`"""
        return DiscreteLpVector


class DiscreteLpVector(DiscretizedSpaceVector):

    """Representation of a `DiscreteLp` element."""

    def asarray(self, out=None):
        """Extract the data of this array as a numpy array.

        Parameters
        ----------
        out : `numpy.ndarray`, optional
            Array in which the result should be written in-place.
            Has to be contiguous and of the correct dtype and
            shape.
        """
        if out is None:
            return super().asarray().reshape(self.shape,
                                             order=self.space.order)
        else:
            if out.shape not in (self.space.shape, (self.space.size,)):
                raise ValueError('output array has shape {}, expected '
                                 '{} or ({},)'
                                 ''.format(out.shape, self.space.shape,
                                           self.space.size))
            out_r = out.reshape(self.space.shape,
                                order=self.space.order)
            if out_r.flags.c_contiguous:
                out_order = 'C'
            elif out_r.flags.f_contiguous:
                out_order = 'F'
            else:
                raise ValueError('output array not contiguous')

            if out_order != self.space.order:
                raise ValueError('output array has ordering {!r}, '
                                 'expected {!r}'
                                 ''.format(self.space.order, out_order))

            super().asarray(out=out.ravel(order=self.space.order))
            return out

    @property
    def ndim(self):
        """Number of dimensions."""
        return self.space.ndim

    @property
    def shape(self):
        """Multi-dimensional shape of this discrete function."""
        # override shape
        return self.space.shape

    @property
    def cell_sides(self):
        """Side lengths of a cell in an underlying *uniform* partition."""
        return self.space.cell_sides

    @property
    def cell_volume(self):
        """Cell volume of an underlying regular grid."""
        return self.space.cell_volume

    @property
    def order(self):
        """Axis ordering for array flattening."""
        return self.space.order

    @property
    def real(self):
        """Real part of this element."""
        rspace = self.space.astype(self.space._real_dtype)
        return rspace.element(self.asarray().real)

    @real.setter
    def real(self, newreal):
        """Set the real part of this element to ``newreal``."""
        # rspace = self.space.astype(self.space._real_dtype)
        # newreal = rspace.element(newreal)
        # self.ntuple.real = newreal.ntuple
        newreal_flat = np.asarray(newreal, order=self.space.order).reshape(
            -1, order=self.space.order)
        self.ntuple.real = newreal_flat

    @property
    def imag(self):
        """Imaginary part of this element."""
        rspace = self.space.astype(self.space._real_dtype)
        return rspace.element(self.asarray().imag)

    @imag.setter
    def imag(self, newimag):
        """Set the imaginary part of this element to ``newimag``."""
        newimag_flat = np.asarray(newimag, order=self.space.order).reshape(
            -1, order=self.space.order)
        self.ntuple.imag = newimag_flat

    def conj(self, out=None):
        """The complex conjugate of this element.

        Parameters
        ----------
        out : `DiscreteLpVector`, optional
            Element to which the complex conjugate is written.
            Must be an element of this vector's space.

        Returns
        -------
        out : `DiscreteLpVector`
            The complex conjugate vector. If ``out`` is provided,
            the returned object is a reference to it.

        Examples
        --------
        >>> discr = uniform_discr(0, 1, 4, dtype='complex')
        >>> x = discr.element([5+1j, 3, 2-2j, 1j])
        >>> y = x.conj(); print(y)
        [(5-1j), (3-0j), (2+2j), -1j]

        The out parameter allows you to avoid a copy:

        >>> z = discr.element()
        >>> z_out = x.conj(out=z); print(z)
        [(5-1j), (3-0j), (2+2j), -1j]
        >>> z_out is z
        True

        It can also be used for in-place conjugation:

        >>> x_out = x.conj(out=x); print(x)
        [(5-1j), (3-0j), (2+2j), -1j]
        >>> x_out is x
        True
        """
        if out is None:
            return self.space.element(self.ntuple.conj())
        else:
            self.ntuple.conj(out=out.ntuple)
            return out

    def __setitem__(self, indices, values):
        """Set values of this vector.

        Parameters
        ----------
        indices : `int` or `slice`
            The position(s) that should be set
        values : scalar, `array-like` or `NtuplesVector`
            The value(s) that are to be assigned.
            If ``indices`` is an `int`, ``values`` must be a single
            value.
            If ``indices`` is a `slice`, ``values`` must be
            broadcastable to the size of the slice (same size,
            shape ``(1,)`` or single value).
            For ``indices=slice(None)``, i.e. in the call
            ``vec[:] = values``, a multi-dimensional array of correct
            shape is allowed as ``values``.
        """
        if values in self.space:
            # For DiscretizedSetVector of the same type, use ntuple directly
            self.ntuple[indices] = values.ntuple
        else:
            # Other sequence types are piped through a Numpy array. Equivalent
            # views are optimized for in Numpy.
            if indices == slice(None):
                values = np.atleast_1d(values)
                if (values.ndim > 1 and
                        values.shape != self.space.shape):
                    raise ValueError('shape {} of value array {} not equal '
                                     'to sampling grid shape {}'
                                     ''.format(values.shape, values,
                                               self.space.shape))
                values = values.ravel(order=self.space.order)

            super().__setitem__(indices, values)

    @property
    def ufunc(self):
        """`DiscreteLpUFuncs`, access to numpy style ufuncs.

        Examples
        --------
        >>> X = uniform_discr(0, 1, 2)
        >>> x = X.element([1, -2])
        >>> x.ufunc.absolute()
        uniform_discr(0.0, 1.0, 2).element([1.0, 2.0])

        These functions can also be used with broadcasting

        >>> x.ufunc.add(3)
        uniform_discr(0.0, 1.0, 2).element([4.0, 1.0])

        and non-space elements

        >>> x.ufunc.subtract([3, 3])
        uniform_discr(0.0, 1.0, 2).element([-2.0, -5.0])

        There is also support for various reductions (sum, prod, min, max)

        >>> x.ufunc.sum()
        -1.0

        Also supports out parameter

        >>> y = X.element([3, 4])
        >>> out = X.element()
        >>> result = x.ufunc.add(y, out=out)
        >>> result
        uniform_discr(0.0, 1.0, 2).element([4.0, 2.0])
        >>> result is out
        True

        Notes
        -----
        These are optimized to use the underlying ntuple space and incur no
        overhead unless these do.
        """
        return DiscreteLpUFuncs(self)

    def show(self, title=None, method='', coords=None, indices=None,
             show=False, fig=None, **kwargs):
        """Display the function graphically.

        Parameters
        ----------
        title : `str`, optional
            Set the title of the figure

        method : `str`, optional
            1d methods:

            'plot' : graph plot

            'scatter' : scattered 2d points
            (2nd axis <-> value)

            2d methods:

            'imshow' : image plot with coloring according to value,
            including a colorbar.

            'scatter' : cloud of scattered 3d points
            (3rd axis <-> value)

        coords : array-like, optional
            Display a slice of the array instead of the full array. The values
            are shown accordinging to the given values, `None` represent all
            values along that dimension. For example [None, None, 0.5] shows
            all values in the first two dimensions, with third coordinate equal
            to 0.5.
            This option is mutually exclusive to indices.

        indices : index expression, optional
            Display a slice of the array instead of the full array. The
            index expression is most easily created with the `numpy.s_`
            constructor, i.e. supply ``np.s_[:, 1, :]`` to display the
            first slice along the second axis.
            For data with 3 or more dimensions, the 2d slice in the first
            two axes at the "middle" along the remaining axes is shown
            (semantically ``[:, :, shape[2:] // 2]``).
            This option is mutually exclusive to coords.

        show : `bool`, optional
            If the plot should be showed now or deferred until later.

        fig : `matplotlib.figure.Figure`
            The figure to show in. Expected to be of same "style", as
            the figure given by this function. The most common use case
            is that ``fig`` is the return value from an earlier call to
            this function.

        kwargs : {'figsize', 'saveto', 'clim', ...}
            Extra keyword arguments passed on to display method
            See the Matplotlib functions for documentation of extra
            options.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure. It is also shown to the user.

        See Also
        --------
        odl.util.graphics.show_discrete_data : Underlying implementation
        """

        from odl.util.graphics import show_discrete_data

        if coords is not None:
            if indices is not None:
                raise ValueError('cannot provide both coords and indices')

            interv = self.space.uspace.domain
            shape = self.shape
            indices = []
            for begin, end, n, coord in zip(interv.begin, interv.end,
                                            shape, coords):
                if coord is None:
                    indices += [slice(None)]
                else:
                    if not begin <= coord <= end:
                        raise ValueError('coord {} not in range {}'
                                         ''.format(coord, [begin, end]))

                    if begin == end:
                        indices += [0]
                    else:
                        normalized_pos = (float(coord) - begin) / (end - begin)
                        indices += [int(n * normalized_pos)]

        # Default to showing x-y slice "in the middle"
        if indices is None and self.ndim >= 3:
            indices = [np.s_[:]] * 2
            indices += [n // 2 for n in self.space.shape[2:]]

        if isinstance(indices, (Integral, slice)):
            indices = (indices,)
        elif indices is None or indices == Ellipsis:
            indices = (np.s_[:],) * self.ndim
        else:
            indices = tuple(indices)

        if Ellipsis in indices:
            # Replace Ellipsis with the correct number of [:] expressions
            pos = indices.index(Ellipsis)
            indices = (indices[:pos] +
                       (np.s_[:], ) * (self.ndim - len(indices) + 1) +
                       indices[pos + 1:])

        if len(indices) < self.ndim:
            raise ValueError('too few axes ({} < {})'.format(len(indices),
                                                             self.ndim))
        if len(indices) > self.ndim:
            raise ValueError('too many axes ({} > {})'.format(len(indices),
                                                              self.ndim))

        if self.ndim <= 3:
            axis_labels = ['x', 'y', 'z']
        else:
            axis_labels = ['x{}'.format(axis) for axis in range(self.ndim)]
        squeezed_axes = [axis for axis in range(self.ndim)
                         if not isinstance(indices[axis], Integral)]
        axis_labels = [axis_labels[axis] for axis in squeezed_axes]

        # Squeeze grid and values according to the index expression
        grid = self.space.grid[indices].squeeze()
        values = self.asarray()[indices].squeeze()

        return show_discrete_data(values, grid, title=title, method=method,
                                  show=show, fig=fig, axis_labels=axis_labels,
                                  **kwargs)


def uniform_discr_frompartition(partition, exponent=2.0, interp='nearest',
                                impl='numpy', **kwargs):
    """Discretize an Lp function space given a uniform partition.

    Parameters
    ----------
    partition : `RectPartition`
        Regular (uniform) partition to be used for discretization
    exponent : positive `float`, optional
        The parameter ``p`` in ``L^p``. If the exponent is not
        equal to the default 2.0, the space has no inner product.
    interp : `str` or `sequence` of `str`, optional
        Interpolation type to be used for discretization.
        A sequence is interpreted as interpolation scheme per axis.

            'nearest' : use nearest-neighbor interpolation

            'linear' : use linear interpolation

    impl : {'numpy', 'cuda'}, optional
        Implementation of the data storage arrays

    Other Parameters
    ----------------
    order : {'C', 'F'}, optional
        Axis ordering in the data storage. Default: 'C'
    dtype : dtype
        Data type for the discretized space

            Default for 'numpy': 'float64' / 'complex128'

            Default for 'cuda': 'float32'

    weighting : {'const', 'none'}, optional
        Weighting of the discretized space functions.

            'const' : weight is a constant, the cell volume (default)

            'none' : no weighting

    Returns
    -------
    discr : `DiscreteLp`
        The uniformly discretized function space

    Examples
    --------
    >>> from odl import uniform_partition
    >>> part = uniform_partition(0, 1, 10)
    >>> uniform_discr_frompartition(part)
    uniform_discr(0.0, 1.0, 10)

    See also
    --------
    uniform_discr : implicit uniform Lp discretization
    uniform_discr_fromspace : uniform Lp discretization from an existing
        function space
    odl.discr.partition.uniform_partition :
        partition of the function domain
    """
    if not isinstance(partition, RectPartition):
        raise TypeError('`partition` {!r} is not a `RectPartition` instance'
                        ''.format(partition))
    if not partition.is_uniform:
        raise ValueError('`partition` is not uniform')

    impl, impl_in = str(impl).lower(), impl
    if impl == 'numpy':
        dtype = np.dtype(kwargs.pop('dtype', 'float64'))
    elif impl == 'cuda':
        if not CUDA_AVAILABLE:
            raise ValueError('cUDA not available')
        dtype = np.dtype(kwargs.pop('dtype', 'float32'))
    else:
        raise ValueError("`impl` '{}' not understood".format(impl_in))

    fspace = FunctionSpace(partition.set, out_dtype=dtype)
    ds_type = dspace_type(fspace, impl, dtype)

    order = kwargs.pop('order', 'C')

    weighting = kwargs.pop('weighting', 'const')
    weighting, weighting_in = str(weighting).lower(), weighting
    if weighting == 'none' or float(exponent) == float('inf'):
        weight = None
    elif weighting == 'const':
        weight = partition.cell_volume
    else:
        raise ValueError("`weighting` '{}' not understood"
                         "".format(weighting_in))

    if dtype is not None:
        dspace = ds_type(partition.size, dtype=dtype, weight=weight,
                         exponent=exponent)
    else:
        dspace = ds_type(partition.size, weight=weight, exponent=exponent)

    return DiscreteLp(fspace, partition, dspace, exponent, interp, order=order)


def uniform_discr_fromspace(fspace, nsamples, exponent=2.0, interp='nearest',
                            impl='numpy', **kwargs):
    """Discretize an Lp function space by uniform partition.

    Parameters
    ----------
    fspace : `FunctionSpace`
        Continuous function space. Its domain must be an
        `IntervalProd` instance.
    nsamples : `int` or `tuple` of `int`
        Number of samples per axis. For dimension >= 2, a tuple is
        required.
    exponent : positive `float`, optional
        The parameter ``p`` in ``L^p``. If the exponent is not
        equal to the default 2.0, the space has no inner product.
    interp : `str` or `sequence` of `str`, optional
        Interpolation type to be used for discretization.
        A sequence is interpreted as interpolation scheme per axis.

            'nearest' : use nearest-neighbor interpolation

            'linear' : use linear interpolation

    impl : {'numpy', 'cuda'}, optional
        Implementation of the data storage arrays

    Other Parameters
    ----------------
    nodes_on_bdry : `bool` or boolean `array-like`, optional
        If `True`, place the outermost grid points at the boundary. For
        `False`, they are shifted by half a cell size to the 'inner'.
        If an array-like is given, it must have shape ``(ndim, 2)``,
        where ``ndim`` is the number of dimensions. It defines per axis
        whether the leftmost (first column) and rightmost (second column)
        nodes node lie on the boundary.
        Default: `False`
    order : {'C', 'F'}, optional
        Axis ordering in the data storage. Default: 'C'
    dtype : dtype, optional
        Data type for the discretized space. If not specified, the
        `FunctionSpace.out_dtype` of ``fspace`` is used.
    weighting : {'const', 'none'}, optional
        Weighting of the discretized space functions.

            'const' : weight is a constant, the cell volume (default)

            'none' : no weighting

    Returns
    -------
    discr : `DiscreteLp`
        The uniformly discretized function space

    Examples
    --------
    >>> from odl import Interval, FunctionSpace
    >>> intv = Interval(0, 1)
    >>> space = FunctionSpace(intv)
    >>> uniform_discr_fromspace(space, 10)
    uniform_discr(0.0, 1.0, 10)

    See also
    --------
    uniform_discr : implicit uniform Lp discretization
    uniform_discr_frompartition : uniform Lp discretization using a given
        uniform partition of a function domain
    uniform_discr_fromintv : uniform discretization from an existing
        interval product
    odl.discr.partition.uniform_partition :
        partition of the function domain
    """
    if not isinstance(fspace, FunctionSpace):
        raise TypeError('`fspace` {!r} is not a `FunctionSpace` instance'
                        ''.format(fspace))
    if not isinstance(fspace.domain, IntervalProd):
        raise TypeError('domain {!r} of the function space is not an '
                        '`IntervalProd` instance'.format(fspace.domain))

    impl, impl_in = str(impl).lower(), impl
    dtype = kwargs.pop('dtype', None)

    # Set data type. If given check consistency with fspace's field and
    # out_dtype. If not given, take the latter.
    if dtype is None:
        dtype = fspace.out_dtype
    else:
        dtype, dtype_in = np.dtype(dtype), dtype
        if not np.can_cast(fspace.out_dtype, dtype, casting='safe'):
            raise ValueError('cannot safely cast from output data {} type of '
                             'the function space to given data type {}'
                             ''.format(fspace.out, dtype_in))

    if fspace.field == RealNumbers() and not is_real_dtype(dtype):
        raise ValueError('cannot discretize real space {} with '
                         'non-real data type {}'
                         ''.format(fspace, dtype))
    elif (fspace.field == ComplexNumbers() and
          not is_complex_floating_dtype(dtype)):
        raise ValueError('cannot discretize complex space {} with '
                         'non-complex-floating data type {}'
                         ''.format(fspace, dtype))

    nodes_on_bdry = kwargs.pop('nodes_on_bdry', False)
    partition = uniform_partition_fromintv(fspace.domain, nsamples,
                                           nodes_on_bdry)

    return uniform_discr_frompartition(partition, exponent, interp, impl_in,
                                       dtype=dtype, **kwargs)


def uniform_discr_fromintv(interval, nsamples, exponent=2.0, interp='nearest',
                           impl='numpy', **kwargs):
    """Discretize an Lp function space by uniform sampling.

    Parameters
    ----------
    domain : `IntervalProd`
        The domain of the uniformly discretized space.
    nsamples : `float` or `tuple` of `float`
        Minimum corner of the result.
    nsamples : `int` or `tuple` of `int`
        Number of samples per axis. For dimension >= 2, a tuple is
        required.
    exponent : positive `float`, optional
        The parameter :math:`p` in :math:`L^p`. If the exponent is not
        equal to the default 2.0, the space has no inner product.
    interp : `str` or `sequence` of `str`, optional
        Interpolation type to be used for discretization.
        A sequence is interpreted as interpolation scheme per axis.

            'nearest' : use nearest-neighbor interpolation

            'linear' : use linear interpolation

    impl : {'numpy', 'cuda'}, optional
        Implementation of the data storage arrays
    nodes_on_bdry : `bool` or `sequence`, optional
        If a sequence is provided, it determines per axis whether to
        place the last grid point on the boundary (True) or shift it
        by half a cell size into the interior (False). In each axis,
        an entry may consist in a single `bool` or a 2-tuple of
        `bool`. In the latter case, the first tuple entry decides for
        the left, the second for the right boundary. The length of the
        sequence must be ``array.ndim``.

        A single boolean is interpreted as a global choice for all
        boundaries.
        Default: `False`

    dtype : dtype, optional
        Data type for the discretized space

            Default for 'numpy': 'float64' / 'complex128'

            Default for 'cuda': 'float32'

    order : {'C', 'F'}, optional
        Ordering of the axes in the data storage. 'C' means the
        first axis varies slowest, the last axis fastest;
        vice versa for 'F'.
        Default: 'C'
    weighting : {'const', 'none'}, optional
        Weighting of the discretized space functions.

            'const' : weight is a constant, the cell volume (default)

            'none' : no weighting

    Returns
    -------
    discr : `DiscreteLp`
        The uniformly discretized function space

    Examples
    --------
    >>> from odl import Interval
    >>> intv = Interval(0, 1)
    >>> uniform_discr_fromintv(intv, 10)
    uniform_discr(0.0, 1.0, 10)

    See also
    --------
    uniform_discr : implicit uniform Lp discretization
    uniform_discr_frompartition : uniform Lp discretization using a given
        uniform partition of a function domain
    uniform_discr_fromspace : uniform discretization from an existing
        function space
    """
    # Select field by dtype
    dtype = kwargs.pop('dtype', None)
    if dtype is None:
        if str(impl).lower() == 'cuda':
            dtype = np.dtype('float32')
        else:
            dtype = np.dtype('float64')

    fspace = FunctionSpace(interval, out_dtype=dtype)
    return uniform_discr_fromspace(fspace, nsamples, exponent, interp, impl,
                                   **kwargs)


def uniform_discr(min_corner, max_corner, nsamples,
                  exponent=2.0, interp='nearest', impl='numpy', **kwargs):
    """Discretize an Lp function space by uniform sampling.

    Parameters
    ----------
    min_corner : `float` or `tuple` of `float`
        Minimum corner of the result.
    max_corner : `float` or `tuple` of `float`
        Minimum corner of the result.
    nsamples : `int` or `tuple` of `int`
        Number of samples per axis. For dimension >= 2, a tuple is
        required.
    exponent : positive `float`, optional
        The parameter :math:`p` in :math:`L^p`. If the exponent is not
        equal to the default 2.0, the space has no inner product.
    interp : `str` or `sequence` of `str`, optional
        Interpolation type to be used for discretization.
        A sequence is interpreted as interpolation scheme per axis.

            'nearest' : use nearest-neighbor interpolation

            'linear' : use linear interpolation

    impl : {'numpy', 'cuda'}, optional
        Implementation of the data storage arrays
    nodes_on_bdry : `bool` or `sequence`, optional
        If a sequence is provided, it determines per axis whether to
        place the last grid point on the boundary (True) or shift it
        by half a cell size into the interior (False). In each axis,
        an entry may consist in a single `bool` or a 2-tuple of
        `bool`. In the latter case, the first tuple entry decides for
        the left, the second for the right boundary. The length of the
        sequence must be ``array.ndim``.

        A single boolean is interpreted as a global choice for all
        boundaries.
        Default: `False`

    dtype : dtype, optional
        Data type for the discretized space

            Default for 'numpy': 'float64' / 'complex128'

            Default for 'cuda': 'float32'

    order : {'C', 'F'}, optional
        Ordering of the axes in the data storage. 'C' means the
        first axis varies slowest, the last axis fastest;
        vice versa for 'F'.
        Default: 'C'
    weighting : {'const', 'none'}, optional
        Weighting of the discretized space functions.

            'const' : weight is a constant, the cell volume (default)

            'none' : no weighting

    Returns
    -------
    discr : `DiscreteLp`
        The uniformly discretized function space

    Examples
    --------
    Create real space:

    >>> uniform_discr([0, 0], [1, 1], [10, 10])
    uniform_discr([0.0, 0.0], [1.0, 1.0], [10, 10])

    Can create complex space by giving a dtype

    >>> uniform_discr([0, 0], [1, 1], [10, 10], dtype='complex')
    uniform_discr([0.0, 0.0], [1.0, 1.0], [10, 10], dtype='complex')

    See also
    --------
    uniform_discr_frompartition : uniform Lp discretization using a given
        uniform partition of a function domain
    uniform_discr_fromspace : uniform discretization from an existing
        function space
    uniform_discr_fromintv : uniform discretization from an existing
        interval product
    """
    interval = IntervalProd(min_corner, max_corner)
    return uniform_discr_fromintv(interval, nsamples, exponent, interp, impl,
                                  **kwargs)


def discr_sequence_space(shape, exponent=2.0, impl='numpy', **kwargs):
    """Return an object mimicing the sequence space ``l^p(R^d)``.

    The returned object is a `DiscreteLp` without sampling and
    interpolation operators. It uses a grid with stride 1 and no
    weighting.

    Parameters
    ----------
    shape : `sequence` of `int`
        Multi-dimensional size of the elements in this space
    exponent : positive `float`, optional
        The parameter ``p`` in ```L^p``. If the exponent is
        not equal to the default 2.0, the space has no inner
        product.
    impl : {'numpy', 'cuda'}
        Implementation of the data storage arrays
    dtype : dtype
        Data type for the discretized space

            Default for 'numpy': 'float64'

            Default for 'cuda': 'float32'

    order : {'C', 'F'}, optional
        Ordering of the axes in the data storage. 'C' means the
        first axis varies slowest, the last axis fastest;
        vice versa for 'F'.
        Default: 'C'

    Returns
    -------
    seqspc : `DiscreteLp`
        The sequence-space-like discrete Lp

    Examples
    --------
    >>> seq_spc = discr_sequence_space((3, 3))
    >>> seq_spc.one().norm() == 3.0
    True
    >>> seq_spc = discr_sequence_space((3, 3), exponent=1)
    >>> seq_spc.one().norm() == 9.0
    True
    """
    kwargs.pop('weighting', None)
    kwargs.pop('nodes_on_bdry', None)
    shape = np.atleast_1d(shape)
    return uniform_discr([0] * len(shape), shape - 1, shape, impl=impl,
                         exponent=exponent, nodes_on_bdry=True,
                         weighting='none', **kwargs)


def uniform_discr_fromdiscr(discr, min_corner=None, max_corner=None,
                            nsamples=None, cell_sides=None, exponent=2.0,
                            interp='nearest', impl='numpy', **kwargs):
    """Return a discretization based on an existing one.

    The parameters that are explicitly given are used to create the
    new discretization, and the missing information is taken from
    the template space. See Notes for the exact functionality.

    Parameters
    ----------
    discr : `DiscreteLp`
        Uniformly discretized space used as a template.
    min_corner : float or sequence of float
        Minimum corner of the resulting spatial domain.
    max_corner : float or sequence of float
        Maximum corner of the resulting spatial domain.
    nsamples : int or sequence of int
        Number of samples per axis.
    exponent : positive float, optional
        The parameter :math:`p` in :math:`L^p`. If the exponent is not
        equal to the default 2.0, the space has no inner product.
    interp : str or sequence of str, optional
        Interpolation type to be used for discretization.
        A sequence is interpreted as interpolation scheme per axis.

            'nearest' : use nearest-neighbor interpolation

            'linear' : use linear interpolation

    impl : {'numpy', 'cuda'}, optional
        Implementation of the data storage arrays
    nodes_on_bdry : bool or sequence, optional
        Specifies whether to put the outmost grid nodes on the
        boundary of the domain.

        If a sequence is provided, it determines per axis whether to
        place the last grid point on the boundary (True) or shift it
        by half a cell size into the interior (False). In each axis,
        an entry may consist in a single boolean or a 2-tuple of
        bool. In the latter case, the first tuple entry decides for
        the left, the second for the right boundary. The length of the
        sequence must be ``discr.ndim``.

        A single boolean is interpreted as a global choice for all
        boundaries.
        Default: `False`

    dtype : optional
        Data type for the discretized space.

            Default for 'numpy': 'float64' / 'complex128'

            Default for 'cuda': 'float32'

    order : {'C', 'F'}, optional
        Ordering of the axes in the data storage. 'C' means the
        first axis varies slowest, the last axis fastest;
        vice versa for 'F'.
        Default: 'C'
    weighting : {'const', 'none'}, optional
        Weighting of the discretized space functions.

            'const' : weight is a constant, the cell volume (default)

            'none' : no weighting

    Notes
    -----
    The parameters ``min_corner``, ``max_corner``, ``nsamples`` and
    ``cell_sides`` can be combined in the following ways (applies in
    each axis individually):

    **0 arguments:**
        Return a copy of ``discr``

    **1 argument:**
        [min,max]_corner -> keep sampling but translate domain so it
        starts/ends at ``[min,max]_corner``

        nsamples/cell_sides -> keep domain but change sampling.
        See `uniform_partition` for restrictions.

    **2 arguments:**
        min_corner + max_corner -> translate and resample with the same
        number of samples

        [min,max]_corner + nsamples/cell_sides -> translate and resample

        nsamples + cell_sides -> error due to ambiguity (keep
        ``min_corner`` or ``max_corner``?)

    **3+ arguments:**
        The underlying partition is uniquely determined by the new
        parameters. See `uniform_partition`.

    See also
    --------
    uniform_partition : underlying domain partitioning scheme

    Examples
    --------
    >>> discr = uniform_discr([0, 0], [1, 2], (10, 5))
    >>> discr.cell_sides
    array([ 0.1,  0.4])

    If no additional argument is given, a copy of ``discr`` is
    returned:

    >>> uniform_discr_fromdiscr(discr) == discr
    True
    >>> uniform_discr_fromdiscr(discr) is discr
    False

    Giving ``min_corner`` or ``max_corner`` results in a
    translation, while for the other two options, the domain
    is kept but re-partitioned:

    >>> uniform_discr_fromdiscr(discr, min_corner=[1, 1])
    uniform_discr([1.0, 1.0], [2.0, 3.0], [10, 5])
    >>> uniform_discr_fromdiscr(discr, max_corner=[0, 0])
    uniform_discr([-1.0, -2.0], [0.0, 0.0], [10, 5])
    >>> uniform_discr_fromdiscr(discr, cell_sides=[1, 1])
    uniform_discr([0.0, 0.0], [1.0, 2.0], [1, 2])
    >>> uniform_discr_fromdiscr(discr, nsamples=[5, 5])
    uniform_discr([0.0, 0.0], [1.0, 2.0], [5, 5])
    >>> uniform_discr_fromdiscr(discr, nsamples=[5, 5]).cell_sides
    array([ 0.2,  0.4])

    The cases with 2 or more additional arguments and the syntax
    for specifying quantities per axis is illustrated in the following:

    # axis 0: translate to match max_corner = 3
    # axis 1: recompute max_corner using the original nsamples with the
    # new min_corner and cell_sides
    >>> new_discr = uniform_discr_fromdiscr(discr, min_corner=[None, 1],
    ...                                     max_corner=[3, None],
    ...                                     cell_sides=[None, 0.25])
    >>> new_discr
    uniform_discr([2.0, 1.0], [3.0, 2.25], [10, 5])
    >>> new_discr.cell_sides
    array([ 0.1 ,  0.25])

    # axis 0: recompute min_corner from old cell_sides and new
    # max_corner and nsamples
    # axis 1: use new min_corner, nsamples and cell_sides only
    >>> new_discr = uniform_discr_fromdiscr(discr, min_corner=[None, 1],
    ...                                     max_corner=[3, None],
    ...                                     nsamples=[5, 5],
    ...                                     cell_sides=[None, 0.25])
    >>> new_discr
    uniform_discr([2.5, 1.0], [3.0, 2.25], [5, 5])
    >>> new_discr.cell_sides
    array([ 0.1 ,  0.25])
    """
    if not isinstance(discr, DiscreteLp):
        raise TypeError('`discr` {!r} is not a DiscreteLp instance'
                        ''.format(discr))
    if not discr.is_uniform:
        raise ValueError('`discr` {} is not uniformly discretized'
                         ''.format(discr))

    # Normalize partition parameters
    min_corner = normalized_scalar_param_list(min_corner, discr.ndim,
                                              param_conv=float, keep_none=True)
    max_corner = normalized_scalar_param_list(max_corner, discr.ndim,
                                              param_conv=float, keep_none=True)
    nsamples = normalized_scalar_param_list(nsamples, discr.ndim,
                                            param_conv=safe_int_conv,
                                            keep_none=True)
    cell_sides = normalized_scalar_param_list(cell_sides, discr.ndim,
                                              param_conv=float, keep_none=True)

    nodes_on_bdry = kwargs.pop('nodes_on_bdry', False)
    nodes_on_bdry = normalized_nodes_on_bdry(nodes_on_bdry, discr.ndim)

    new_beg = []
    new_end = []
    new_shape = []
    new_csides = []
    for i, (b, e, n, s, old_b, old_e, old_n, old_s) in enumerate(
            zip(min_corner, max_corner, nsamples, cell_sides,
                discr.min_corner, discr.max_corner, discr.shape,
                discr.cell_sides)):
        num_params = sum(p is not None for p in (b, e, n, s))

        if num_params == 0:
            new_params = [old_b, old_e, old_n, None]

        elif num_params == 1:
            if b is not None:
                new_params = [b, old_e + (b - old_b), old_n, None]
            elif e is not None:
                new_params = [old_b + (e - old_e), e, old_n, None]
            elif n is not None:
                new_params = [old_b, old_e, n, None]
            else:
                new_params = [old_b, old_e, None, s]

        elif num_params == 2:
            if b is not None and e is not None:
                new_params = [b, e, old_n, None]
            elif b is not None and n is not None:
                new_params = [b, None, n, old_s]
            elif b is not None and s is not None:
                new_params = [b, None, old_n, s]
            elif e is not None and n is not None:
                new_params = [None, e, n, old_s]
            elif e is not None and s is not None:
                new_params = [None, e, old_n, s]
            else:
                raise ValueError('in axis {}: cannot use `nsamples` and '
                                 '`cell_size` only due to ambiguous values '
                                 'for begin and end.'.format(i))

        else:
            new_params = [b, e, n, s]

        new_beg.append(new_params[0])
        new_end.append(new_params[1])
        new_shape.append(new_params[2])
        new_csides.append(new_params[3])

    new_part = uniform_partition(begin=new_beg, end=new_end,
                                 num_nodes=new_shape,
                                 cell_sides=new_csides,
                                 nodes_on_bdry=nodes_on_bdry)

    # TODO: use property when it becomes available
    impl = 'cuda' if isinstance(discr.dspace, CudaFn) else 'numpy'
    return uniform_discr_frompartition(new_part, exponent=discr.exponent,
                                       interp=discr.interp, impl=impl)


def _scaling_func_list(bdry_fracs, exponent=1.0):
    """Return a list of lists of scaling functions for the boundary."""
    def scaling(factor):
        def scaling_func(x):
            return x * factor
        return scaling_func

    func_list = []
    for frac_l, frac_r in bdry_fracs:
        func_list_entry = []
        if np.isclose(frac_l, 1.0):
            func_list_entry.append(None)
        else:
            func_list_entry.append(scaling(frac_l ** (1 / exponent)))

        if np.isclose(frac_r, 1.0):
            func_list_entry.append(None)
        else:
            func_list_entry.append(scaling(frac_r ** (1 / exponent)))

        func_list.append(func_list_entry)
    return func_list


if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
