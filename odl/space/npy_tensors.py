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

"""CPU implementations of tensor spaces."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
from future.utils import native
standard_library.install_aliases()
from builtins import super

import ctypes
from functools import partial
from math import sqrt
from numbers import Integral
import numpy as np
import scipy.linalg as linalg
from scipy.sparse.base import isspmatrix

from odl.operator.operator import Operator
from odl.set.sets import RealNumbers, ComplexNumbers
from odl.space.base_ntuples import (
    NtuplesBase, NtuplesBaseVector, FnBase, FnBaseVector)
from odl.space.base_tensors import (
    TensorSetBase, PreBaseTensor, TensorSpaceBase, BaseTensor)
from odl.space.weighting import (
    WeightingBase, MatrixWeightingBase, VectorWeightingBase,
    ConstWeightingBase, NoWeightingBase,
    CustomInnerProductBase, CustomNormBase, CustomDistBase)
from odl.util.ufuncs import NtuplesUFuncs
from odl.util.utility import (
    dtype_repr, is_real_dtype, is_real_floating_dtype,
    is_complex_floating_dtype)


__all__ = ('NumpyTensorSet', 'NumpyPreTensor',
           'NumpyTensorSpace', 'NumpyTensor')


_BLAS_DTYPES = (np.dtype('float32'), np.dtype('float64'),
                np.dtype('complex64'), np.dtype('complex128'))


class NumpyTensorSet(TensorSetBase):

    """The set of tensors of arbitrary type."""

    def element(self, inp=None, data_ptr=None):
        """Create a new element.

        Parameters
        ----------
        inp : `array-like`, optional
            Input to initialize the new element with.

            If ``inp`` is `None`, an empty element is created with no
            guarantee of its state (memory allocation only).

            If ``inp`` is a `numpy.ndarray` of the same shape and data
            type as this space, the array is wrapped, not copied.
            Other array-like objects are copied.

        Returns
        -------
        element : `NumpyPreTensor`
            The new element created (from ``inp``).

        Notes
        -----
        This method preserves "array views" of correct size and type,
        see the examples below.

        Examples
        --------
        TODO: rewrite
        >>> strings3 = Ntuples(3, dtype='U1')  # 1-char strings
        >>> x = strings3.element(['w', 'b', 'w'])
        >>> print(x)

        >>> x.space


        Construction from data pointer:

        >>> int3 = Ntuples(3, dtype='int')
        >>> x = int3.element([1, 2, 3])
        >>> y = int3.element(data_ptr=x.data_ptr)
        >>> print(y)

        >>> y[0] = 5
        >>> print(x)

        """
        # TODO: exactly the same code would work for NumpyNtuples, factor out!
        if inp is None:
            if data_ptr is None:
                arr = np.empty(self.shape, dtype=self.dtype, order=self.order)
                return self.element_type(self, arr)
            else:
                ctype_array_def = ctypes.c_byte * self.nbytes
                as_ctype_array = ctype_array_def.from_address(data_ptr)
                as_numpy_array = np.ctypeslib.as_array(as_ctype_array)
                arr = as_numpy_array.view(dtype=self.dtype).reshape(self.shape)
                return self.element_type(self, arr)
        else:
            if data_ptr is None:
                if inp in self:
                    return inp
                else:
                    arr = np.array(inp, copy=False, dtype=self.dtype, ndmin=1)
                    if arr.shape != self.shape:
                        raise ValueError('expected input shape {}, got {}'
                                         ''.format((self.size,), arr.shape))

                    return self.element_type(self, arr)
            else:
                raise ValueError('cannot provide both `inp` and `data_ptr`')

    def zero(self):
        """Return a tensor of all zeros.

        Examples
        --------
        >>> c3 = Cn(3)
        >>> x = c3.zero()
        >>> x

        """
        return self.element(np.zeros(self.shape, dtype=self.dtype,
                                     order=self.order))

    def one(self):
        """Return a tensor of all ones.

        Examples
        --------
        >>> c3 = Cn(3)
        >>> x = c3.one()
        >>> x

        """
        return self.element(np.ones(self.shape, dtype=self.dtype,
                                    order=self.order))

    @property
    def element_type(self):
        """Type of elements in this space."""
        return NumpyPreTensor


class NumpyPreTensor(BaseTensor):

    """Representation of a `NumpyTensorSet` element."""

    def __init__(self, space, data):
        """Initialize a new instance."""
        if not isinstance(space, NumpyTensorSet):
            raise TypeError('`space` must be a `NumpyTensorSet` '
                            'instance, got {!r}'.format(space))

        if not isinstance(data, np.ndarray):
            raise TypeError('`data` {!r} not a `numpy.ndarray` instance'
                            ''.format(data))

        if data.dtype != space.dtype:
            raise TypeError('`data` {!r} not of dtype {!r}'
                            ''.format(data, space.dtype))
        self._data = data

        BaseTensor.__init__(self, space)

    @property
    def data(self):
        """The raw `numpy.ndarray` representing the data."""
        return self._data

    def asarray(self, out=None):
        """Extract the data of this array as a ``numpy.ndarray``.

        Parameters
        ----------
        out : `numpy.ndarray`, optional
            Array in which the result should be written in-place.
            Has to be contiguous and of the correct dtype.

        Returns
        -------
        asarray : `numpy.ndarray`
            Numpy array with the same data type as ``self``. If
            ``out`` was given, the returned object is a reference
            to it.

        Examples
        --------
        >>> vec = Ntuples(3, 'float').element([1, 2, 3])
        >>> vec.asarray()

        >>> vec.asarray(start=1, stop=3)


        Using the out parameter:

        >>> out = np.empty((3,), dtype='float')
        >>> result = vec.asarray(out=out)
        >>> out

        >>> result is out

        """
        if out is None:
            return self.data
        else:
            out[:] = self.data
            return out

    @property
    def data_ptr(self):
        """A raw pointer to the data container.

        Examples
        --------
        >>> import ctypes
        >>> vec = Ntuples(3, 'int32').element([1, 2, 3])
        >>> arr_type = ctypes.c_int32 * 3
        >>> buffer = arr_type.from_address(vec.data_ptr)
        >>> arr = np.frombuffer(buffer, dtype='int32')
        >>> print(arr)


        In-place modification via pointer:

        >>> arr[0] = 5
        >>> print(vec)

        """
        return self.data.ctypes.data

    def __eq__(self, other):
        """Return ``self == other``.

        Returns
        -------
        equals : `bool`
            `True` if all entries of ``other`` are equal to this
            the entries of ``self``, `False` otherwise.

        Notes
        -----
        Space membership is not checked, hence vectors from
        different spaces can be equal.

        Examples
        --------
        >>> vec1 = Ntuples(3, int).element([1, 2, 3])
        >>> vec2 = Ntuples(3, int).element([-1, 2, 0])
        >>> vec1 == vec2

        >>> vec2 = Ntuples(3, int).element([1, 2, 3])
        >>> vec1 == vec2


        Space membership matters:

        >>> vec2 = Ntuples(3, float).element([1, 2, 3])
        >>> vec1 == vec2 or vec2 == vec1

        """
        if other is self:
            return True
        elif other not in self.space:
            return False
        else:
            return np.array_equal(self.data, other.data)

    def copy(self):
        """Create an identical (deep) copy of this vector.

        Parameters
        ----------
        None

        Returns
        -------
        copy : `NumpyPreTensor`
            The deep copy

        Examples
        --------
        >>> vec1 = Ntuples(3, 'int').element([1, 2, 3])
        >>> vec2 = vec1.copy()
        >>> vec2

        >>> vec1 == vec2

        >>> vec1 is vec2

        """
        return self.space.element(self.data.copy())

    __copy__ = copy

    def __getitem__(self, indices):
        """Return ``self[indices]``.

        Parameters
        ----------
        indices : index expression
            Integer, slice or sequence of these, defining the positions
            of the data array which should be accessed.

        Returns
        -------
        values : `NumpyTensorSet.dtype` or `NumpyPreTensor`
            The value(s) at the given indices. Note that the returned
            object is a writable view into the original tensor.

        Examples
        --------
        >>> str_3 = Ntuples(3, dtype='U6')  # 6-char unicode
        >>> x = str_3.element(['a', 'Hello!', '0'])
        >>> print(x[0])

        >>> print(x[1:3])

        >>> x[1:3].space

        """
        arr = self.data[indices]
        if np.isscalar(arr):
            return arr
        else:
            return type(self.space)(arr.shape, dtype=self.dtype,
                                    order=self.order).element(arr)

    def __setitem__(self, indices, values):
        """Implement ``self[indices] = values``.

        Parameters
        ----------
        indices : index expression
            Integer, slice or sequence of these, defining the positions
            of the data array which should be written to.
        values : scalar, array-like or `BaseTensor`
            The value(s) that are to be assigned.

            If ``index`` is an integer, ``value`` must be a scalar.

            If ``index`` is a slice or a sequence of slices, ``value``
            must be broadcastable to the shape of the slice.

        Examples
        --------
        >>> int_3 = Ntuples(3, 'int')
        >>> x = int_3.element([1, 2, 3])
        >>> x[0] = 5
        >>> x


        Assignment from array-like structures or another
        vector:

        >>> y = Ntuples(2, 'short').element([-1, 2])
        >>> x[:2] = y
        >>> x

        >>> x[1:3] = [7, 8]
        >>> x

        >>> x[:] = np.array([0, 0, 0])
        >>> x


        Broadcasting is also supported:

        >>> x[1:3] = -2.
        >>> x


        Array views are preserved:

        >>> y = x[::2]  # view into x
        >>> y[:] = -9
        >>> print(y)

        >>> print(x)


        Be aware of unsafe casts and over-/underflows, there
        will be warnings at maximum.

        >>> x = Ntuples(2, 'int8').element([0, 0])
        >>> maxval = 255  # maximum signed 8-bit unsigned int
        >>> x[0] = maxval + 1
        >>> x

        """
        if isinstance(values, NumpyPreTensor):
            self.data[indices] = values.data
        else:
            self.data[indices] = values

    @property
    def ufunc(self):
        """`NtuplesUFuncs`, access to numpy style ufuncs.

        Examples
        --------
        >>> r2 = Rn(2)
        >>> x = r2.element([1, -2])
        >>> x.ufunc.absolute()
        Rn(2).element([1.0, 2.0])

        These functions can also be used with broadcasting

        >>> x.ufunc.add(3)
        Rn(2).element([4.0, 1.0])

        and non-space elements

        >>> x.ufunc.subtract([3, 3])
        Rn(2).element([-2.0, -5.0])

        There is also support for various reductions (sum, prod, min, max)

        >>> x.ufunc.sum()
        -1.0

        They also support an out parameter

        >>> y = r2.element([3, 4])
        >>> out = r2.element()
        >>> result = x.ufunc.add(y, out=out)
        >>> result
        Rn(2).element([4.0, 2.0])
        >>> result is out
        True

        Notes
        -----
        These are optimized for use with ntuples and incur no overhead.
        """
        # TODO: implement
        raise NotImplementedError


def _blas_is_applicable(*args):
    """Whether BLAS routines can be applied or not.

    BLAS routines are available for single and double precision
    float or complex data only. If the arrays are non-contiguous,
    BLAS methods are usually slower, and array-writing routines do
    not work at all. Hence, only contiguous arrays are allowed.

    Parameters
    ----------
    x1,...,xN : `NumpyPreTensor`
        The tensors to be tested for BLAS conformity.
    """
    return (all(x.dtype == args[0].dtype and
                x.dtype in _BLAS_DTYPES and
                x.data.flags.contiguous
                for x in args))


def _lincomb(a, x1, b, x2, out, dtype):
    """Raw linear combination depending on data type."""

    # Shortcut for small problems
    if x1.size < 100:  # small array optimization
        out.data[:] = a * x1.data + b * x2.data
        return

    # Use blas for larger problems
    def fallback_axpy(x1, x2, n, a):
        """Fallback axpy implementation avoiding copy."""
        if a != 0:
            x2 /= a
            x2 += x1
            x2 *= a
        return x2

    def fallback_scal(a, x, n):
        """Fallback scal implementation."""
        x *= a
        return x

    def fallback_copy(x1, x2, n):
        """Fallback copy implementation."""
        x2[...] = x1[...]
        return x2

    if _blas_is_applicable(x1, x2, out):
        axpy, scal, copy = linalg.blas.get_blas_funcs(
            ['axpy', 'scal', 'copy'], arrays=(x1.data, x2.data, out.data))
    else:
        axpy, scal, copy = (fallback_axpy, fallback_scal, fallback_copy)

    if x1 is x2 and b != 0:
        # x1 is aligned with x2 -> out = (a+b)*x1
        _lincomb(a + b, x1, 0, x1, out, dtype)
    elif out is x1 and out is x2:
        # All the vectors are aligned -> out = (a+b)*out
        scal(a + b, out.data, native(out.size))
    elif out is x1:
        # out is aligned with x1 -> out = a*out + b*x2
        if a != 1:
            scal(a, out.data, native(out.size))
        if b != 0:
            axpy(x2.data, out.data, native(out.size), b)
    elif out is x2:
        # out is aligned with x2 -> out = a*x1 + b*out
        if b != 1:
            scal(b, out.data, native(out.size))
        if a != 0:
            axpy(x1.data, out.data, native(out.size), a)
    else:
        # We have exhausted all alignment options, so x1 != x2 != out
        # We now optimize for various values of a and b
        if b == 0:
            if a == 0:  # Zero assignment -> out = 0
                out.data[:] = 0
            else:  # Scaled copy -> out = a*x1
                copy(x1.data, out.data, native(out.size))
                if a != 1:
                    scal(a, out.data, native(out.size))
        else:
            if a == 0:  # Scaled copy -> out = b*x2
                copy(x2.data, out.data, native(out.size))
                if b != 1:
                    scal(b, out.data, native(out.size))

            elif a == 1:  # No scaling in x1 -> out = x1 + b*x2
                copy(x1.data, out.data, native(out.size))
                axpy(x2.data, out.data, native(out.size), b)
            else:  # Generic case -> out = a*x1 + b*x2
                copy(x2.data, out.data, native(out.size))
                if b != 1:
                    scal(b, out.data, native(out.size))
                axpy(x1.data, out.data, native(out.size), a)


class NumpyTensorSpace(TensorSpaceBase, NumpyTensorSet):

    """The space of tensors of given shape.

    This space implements multi-dimensional arrays whose entries are
    elements of a `Field`, which is usually the real or complex numbers.

    The space elements are represented as instances of the
    `NumpyTensor` class.
    """

    def __init__(self, shape, dtype, order='C', **kwargs):
        """Initialize a new instance.

        Parameters
        ----------
        size : positive `int`
            The number of dimensions of the space
        dtype : `object`
            The data type of the storage array. Can be provided in any
            way the `numpy.dtype` function understands, most notably
            as built-in type, as `numpy.dtype` or as `string`.

            Only scalar data types are allowed.

        order : {'C', 'F'}, optional
            Axis ordering of the data storage.

        exponent : positive `float`, optional
            Exponent of the norm. For values other than 2.0, no
            inner product is defined.

            This option is ignored if ``dist``, ``norm`` or
            ``inner`` is given.

            Default: 2.0

        Other Parameters
        ----------------

        dist : `callable`, optional
            Distance function defining a metric on the space.
            It must accept two `NumpyTensor` arguments and return
            a non-negative real number. See ``Notes`` for
            mathematical requirements.

            By default, ``dist(x, y)`` is calculated as ``norm(x - y)``.
            This creates an intermediate array ``x - y``, which can be
            avoided by choosing ``dist_using_inner=True``.

            This option cannot be combined with ``weight``,
            ``norm`` or ``inner``.

        norm : `callable`, optional
            The norm implementation. It must accept a
            `NumpyTensor` argument, return a non-negative real number.
            See ``Notes`` for mathematical requirements.

            By default, ``norm(x)`` is calculated as ``inner(x, x)``.

            This option cannot be combined with ``weight``,
            ``dist`` or ``inner``.

        inner : `callable`, optional
            The inner product implementation. It must accept two
            `NumpyTensor` arguments and return an element of the field
            of the space (usually real or complex number).
            See ``Notes`` for mathematical requirements.

            This option cannot be combined with ``weight``,
            ``dist`` or ``norm``.

        dist_using_inner : `bool`, optional
            Calculate ``dist`` using the formula

                ``||x - y||^2 = ||x||^2 + ||y||^2 - 2 * Re <x, y>``

            This avoids the creation of new arrays and is thus faster
            for large arrays. On the downside, it will not evaluate to
            exactly zero for equal (but not identical) ``x`` and ``y``.

            This option can only be used if ``exponent`` is 2.0.

            Default: `False`.

        kwargs :
            Further keyword arguments are passed to the weighting
            classes.

        See also
        --------
        NumpyTensorSpaceVectorWeighting
        NumpyTensorSpaceConstWeighting

        Notes
        -----
        - A distance function or metric on a space :math:`\mathcal{X}`
          is a mapping
          :math:`d:\mathcal{X} \\times \mathcal{X} \\to \mathbb{R}`
          satisfying the following conditions for all space elements
          :math:`x, y, z`:

          * :math:`d(x, y) \geq 0`,
          * :math:`d(x, y) = 0 \Leftrightarrow x = y`,
          * :math:`d(x, y) = d(y, x)`,
          * :math:`d(x, y) \\leq d(x, z) + d(z, y)`.

        - A norm on a space :math:`\mathcal{X}` is a mapping
          :math:`\\lVert \cdot \\rVert:\mathcal{X} \\to \mathbb{R}`
          satisfying the following conditions for all
          space elements :math:`x, y`: and scalars :math:`s`:

          * :math:`\\lVert x\\rVert \geq 0`,
          * :math:`\\lVert x\\rVert = 0 \Leftrightarrow x = 0`,
          * :math:`\\lVert sx\\rVert = |s| \cdot \\lVert x \\rVert`,
          * :math:`\\lVert x+y\\rVert \\leq \\lVert x\\rVert +
            \\lVert y\\rVert`.

        - An inner product on a space :math:`\mathcal{X}` over a field
          :math:`\mathbb{F} = \mathbb{R}` or :math:`\mathbb{C}` is a
          mapping
          :math:`\\langle\cdot, \cdot\\rangle: \mathcal{X} \\times
          \mathcal{X} \\to \mathbb{F}`
          satisfying the following conditions for all
          space elements :math:`x, y, z`: and scalars :math:`s`:

          * :math:`\\langle x, y\\rangle =
            \overline{\\langle y, x\\rangle}`,
          * :math:`\\langle sx + y, z\\rangle = s \\langle x, z\\rangle +
            \\langle y, z\\rangle`,
          * :math:`\\langle x, x\\rangle = 0 \Leftrightarrow x = 0`.

        Examples
        --------
        >>> space = Fn(3, 'float')
        >>> space

        >>> space = Fn(3, 'float', weight=[1, 2, 3])
        >>> space

        """
        # TODO: continue here
        Ntuples.__init__(self, size, dtype)
        FnBase.__init__(self, size, dtype)

        dist = kwargs.pop('dist', None)
        norm = kwargs.pop('norm', None)
        inner = kwargs.pop('inner', None)
        weight = kwargs.pop('weight', None)
        exponent = kwargs.pop('exponent', 2.0)
        dist_using_inner = bool(kwargs.pop('dist_using_inner', False))

        # Check validity of option combination (3 or 4 out of 4 must be None)
        if sum(x is None for x in (dist, norm, inner, weight)) < 3:
            raise ValueError('invalid combination of options `weight`, '
                             '`dist`, `norm` and `inner`')

        if any(x is not None for x in (dist, norm, inner)) and exponent != 2.0:
            raise ValueError('`exponent` cannot be used together with '
                             '`dist`, `norm` and `inner`')

        # Set the weighting
        if weight is not None:
            if isinstance(weight, WeightingBase):
                self._weighting = weight
            elif np.isscalar(weight):
                self._weighting = FnConstWeighting(
                    weight, exponent, dist_using_inner=dist_using_inner)
            elif weight is None:
                # Need to wait until dist, norm and inner are handled
                pass
            elif isspmatrix(weight):
                self._weighting = FnMatrixWeighting(
                    weight, exponent, dist_using_inner=dist_using_inner,
                    **kwargs)
            else:  # last possibility: make a matrix
                arr = np.asarray(weight)
                if arr.dtype == object:
                    raise ValueError('invalid weight argument {}'
                                     ''.format(weight))
                if arr.ndim == 1:
                    self._weighting = FnVectorWeighting(
                        arr, exponent, dist_using_inner=dist_using_inner)
                elif arr.ndim == 2:
                    self._weighting = FnMatrixWeighting(
                        arr, exponent, dist_using_inner=dist_using_inner,
                        **kwargs)
                else:
                    raise ValueError('array-like input {} is not 1- or '
                                     '2-dimensional'.format(weight))

        elif dist is not None:
            self._weighting = FnCustomDist(dist)
        elif norm is not None:
            self._weighting = FnCustomNorm(norm)
        elif inner is not None:
            self._weighting = FnCustomInnerProduct(inner)
        else:  # all None -> no weighing
            self._weighting = FnNoWeighting(
                exponent, dist_using_inner=dist_using_inner)

    @property
    def exponent(self):
        """Exponent of the norm and distance."""
        return self._weighting.exponent

    @property
    def weighting(self):
        """This space's weighting scheme."""
        return self._weighting

    @property
    def is_weighted(self):
        """Return `True` if the weighting is not `FnNoWeighting`."""
        return not isinstance(self.weighting, FnNoWeighting)

    @staticmethod
    def default_dtype(field):
        """Return the default of `Fn` data type for a given field.

        Parameters
        ----------
        field : `Field`
            Set of numbers to be represented by a data type.
            Currently supported : `RealNumbers`, `ComplexNumbers`

        Returns
        -------
        dtype : `type`
            Numpy data type specifier. The returned defaults are:

            ``RealNumbers()`` : ``np.dtype('float64')``

            ``ComplexNumbers()`` : ``np.dtype('complex128')``
        """
        if field == RealNumbers():
            return np.dtype('float64')
        elif field == ComplexNumbers():
            return np.dtype('complex128')
        else:
            raise ValueError('no default data type defined for field {}'
                             ''.format(field))

    def _lincomb(self, a, x1, b, x2, out):
        """Linear combination of ``x1`` and ``x2``.

        Calculate ``out = a*x1 + b*x2`` using optimized BLAS
        routines if possible.

        Parameters
        ----------
        a, b : `FnBase.field`
            Scalars to multiply ``x1`` and ``x2`` with
        x1, x2 : `FnVector`
            Summands in the linear combination
        out : `FnVector`
            Vector to which the result is written

        Returns
        -------
        `None`

        Examples
        --------
        >>> c3 = Cn(3)
        >>> x = c3.element([1+1j, 2-1j, 3])
        >>> y = c3.element([4+0j, 5, 6+0.5j])
        >>> out = c3.element()
        >>> c3.lincomb(2j, x, 3-1j, y, out)  # out is returned
        Cn(3).element([(10-2j), (17-1j), (18.5+1.5j)])
        >>> out
        Cn(3).element([(10-2j), (17-1j), (18.5+1.5j)])
        """
        _lincomb(a, x1, b, x2, out, self.dtype)

    def _dist(self, x1, x2):
        """Calculate the distance between two vectors.

        Parameters
        ----------
        x1, x2 : `FnVector`
            Vectors whose mutual distance is calculated

        Returns
        -------
        dist : `float`
            Distance between the vectors

        Examples
        --------
        >>> from numpy.linalg import norm
        >>> c2_2 = Cn(2, dist=lambda x, y: norm(x - y, ord=2))
        >>> x = c2_2.element([3+1j, 4])
        >>> y = c2_2.element([1j, 4-4j])
        >>> c2_2.dist(x, y)
        5.0

        >>> c2_2 = Cn(2, dist=lambda x, y: norm(x - y, ord=1))
        >>> x = c2_2.element([3+1j, 4])
        >>> y = c2_2.element([1j, 4-4j])
        >>> c2_2.dist(x, y)
        7.0
        """
        return self.weighting.dist(x1, x2)

    def _norm(self, x):
        """Calculate the norm of a vector.

        Parameters
        ----------
        x : `FnVector`
            The vector whose norm is calculated

        Returns
        -------
        norm : `float`
            Norm of the vector

        Examples
        --------
        >>> import numpy as np
        >>> c2_2 = Cn(2, norm=np.linalg.norm)  # 2-norm
        >>> x = c2_2.element([3+1j, 1-5j])
        >>> c2_2.norm(x)
        6.0

        >>> from functools import partial
        >>> c2_1 = Cn(2, norm=partial(np.linalg.norm, ord=1))
        >>> x = c2_1.element([3-4j, 12+5j])
        >>> c2_1.norm(x)
        18.0
        """
        return self.weighting.norm(x)

    def _inner(self, x1, x2):
        """Raw inner product of two vectors.

        Parameters
        ----------
        x1, x2 : `FnVector`
            The vectors whose inner product is calculated

        Returns
        -------
        inner : `field` `element`
            Inner product of the vectors

        Examples
        --------
        >>> import numpy as np
        >>> c3 = Cn(2, inner=lambda x, y: np.vdot(y, x))
        >>> x = c3.element([5+1j, -2j])
        >>> y = c3.element([1, 1+1j])
        >>> c3.inner(x, y) == (5+1j)*1 + (-2j)*(1-1j)
        True

        Define a space with custom inner product:

        >>> weights = np.array([1., 2.])
        >>> def weighted_inner(x, y):
        ...     return np.vdot(weights * y.data, x.data)

        >>> c3w = Cn(2, inner=weighted_inner)
        >>> x = c3w.element(x)  # elements must be cast (no copy)
        >>> y = c3w.element(y)
        >>> c3w.inner(x, y) == 1*(5+1j)*1 + 2*(-2j)*(1-1j)
        True
        """
        return self.weighting.inner(x1, x2)

    def _multiply(self, x1, x2, out):
        """The entry-wise product of two vectors, assigned to out.

        Parameters
        ----------
        x1, x2 : `FnVector`
            Factors in the product
        out : `FnVector`
            Vector to which the result is written

        Returns
        -------
        `None`

        Examples
        --------
        >>> c3 = Cn(3)
        >>> x = c3.element([5+1j, 3, 2-2j])
        >>> y = c3.element([1, 2+1j, 3-1j])
        >>> out = c3.element()
        >>> c3.multiply(x, y, out)  # out is returned
        Cn(3).element([(5+1j), (6+3j), (4-8j)])
        >>> out
        Cn(3).element([(5+1j), (6+3j), (4-8j)])
        """
        np.multiply(x1.data, x2.data, out=out.data)

    def _divide(self, x1, x2, out):
        """The entry-wise division of two vectors, assigned to out.

        Parameters
        ----------
        x1, x2 : `FnVector`
            Dividend and divisor in the quotient
        out : `FnVector`
            Vector to which the result is written

        Returns
        -------
        `None`

        Examples
        --------
        >>> r3 = Rn(3)
        >>> x = r3.element([3, 5, 6])
        >>> y = r3.element([1, 2, 2])
        >>> out = r3.element()
        >>> r3.divide(x, y, out)  # out is returned
        Rn(3).element([3.0, 2.5, 3.0])
        >>> out
        Rn(3).element([3.0, 2.5, 3.0])
        """
        np.divide(x1.data, x2.data, out=out.data)

    def __eq__(self, other):
        """Return ``self == other``.

        Returns
        -------
        equals : `bool`
            `True` if other is an instance of this space's type
            with the same
            `NtuplesBase.size` and `NtuplesBase.dtype`,
            and identical distance function, otherwise `False`.

        Examples
        --------
        >>> from numpy.linalg import norm
        >>> def dist(x, y, ord):
        ...     return norm(x - y, ord)

        >>> from functools import partial
        >>> dist2 = partial(dist, ord=2)
        >>> c3 = Cn(3, dist=dist2)
        >>> c3_same = Cn(3, dist=dist2)
        >>> c3  == c3_same
        True

        Different ``dist`` functions result in different spaces - the
        same applies for ``norm`` and ``inner``:

        >>> dist1 = partial(dist, ord=1)
        >>> c3_1 = Cn(3, dist=dist1)
        >>> c3_2 = Cn(3, dist=dist2)
        >>> c3_1 == c3_2
        False

        Be careful with Lambdas - they result in non-identical function
        objects:

        >>> c3_lambda1 = Cn(3, dist=lambda x, y: norm(x-y, ord=1))
        >>> c3_lambda2 = Cn(3, dist=lambda x, y: norm(x-y, ord=1))
        >>> c3_lambda1 == c3_lambda2
        False

        An `Fn` space with the same data type is considered
        equal:

        >>> c3 = Cn(3)
        >>> f3_cdouble = Fn(3, dtype='complex128')
        >>> c3 == f3_cdouble
        True
        """
        if other is self:
            return True

        return (Ntuples.__eq__(self, other) and
                self.weighting == other.weighting)

    def __repr__(self):
        """Return ``repr(self)``."""
        if self.is_rn:
            class_name = 'Rn'
        elif self.is_cn:
            class_name = 'Cn'
        else:
            class_name = self.__class__.__name__

        inner_str = '{}'.format(self.size)
        if self.dtype != self.default_dtype(self.field):
            inner_str += ', {}'.format(dtype_repr(self.dtype))

        weight_str = self.weighting.repr_part
        if weight_str:
            inner_str += ', ' + weight_str
        return '{}({})'.format(class_name, inner_str)

    # Copy these to handle bug in ABCmeta
    zero = Ntuples.zero
    one = Ntuples.one

    @property
    def element_type(self):
        """Return `FnVector`."""
        return FnVector


class FnVector(FnBaseVector, NtuplesVector):

    """Representation of an `Fn` element."""

    def __init__(self, space, data):
        """Initialize a new instance."""
        if not isinstance(space, Fn):
            raise TypeError('{!r} not an `Fn` instance'
                            ''.format(space))

        FnBaseVector.__init__(self, space)
        NtuplesVector.__init__(self, space, data)

    @property
    def real(self):
        """The real part of this vector.

        Returns
        -------
        real : `FnVector` view with dtype
            The real part this vector as a vector in `Rn`

        Examples
        --------
        >>> c3 = Cn(3)
        >>> x = c3.element([5+1j, 3, 2-2j])
        >>> x.real
        Rn(3).element([5.0, 3.0, 2.0])

        The `Rn` vector is really a view, so changes affect
        the original array:

        >>> x.real *= 2
        >>> x
        Cn(3).element([(10+1j), (6+0j), (4-2j)])
        """
        rn = Rn(self.space.size, self.space._real_dtype)
        return rn.element(self.data.real)

    @real.setter
    def real(self, newreal):
        """The setter for the real part.

        This method is invoked by ``vec.real = other``.

        Parameters
        ----------
        newreal : `array-like` or scalar
            The new real part for this vector.

        Examples
        --------
        >>> c3 = Cn(3)
        >>> x = c3.element([5+1j, 3, 2-2j])
        >>> a = Rn(3).element([0, 0, 0])
        >>> x.real = a
        >>> x
        Cn(3).element([1j, 0j, -2j])

        Other array-like types and broadcasting:

        >>> x.real = 1.0
        >>> x
        Cn(3).element([(1+1j), (1+0j), (1-2j)])
        >>> x.real = [0, 2, -1]
        >>> x
        Cn(3).element([1j, (2+0j), (-1-2j)])
        """
        self.real.data[:] = newreal

    @property
    def imag(self):
        """The imaginary part of this vector.

        Returns
        -------
        imag : `FnVector`
            The imaginary part this vector as a vector in
            `Rn`

        Examples
        --------
        >>> c3 = Cn(3)
        >>> x = c3.element([5+1j, 3, 2-2j])
        >>> x.imag
        Rn(3).element([1.0, 0.0, -2.0])

        The `Rn` vector is really a view, so changes affect
        the original array:

        >>> x.imag *= 2
        >>> x
        Cn(3).element([(5+2j), (3+0j), (2-4j)])
        """
        rn = Rn(self.space.size, self.space._real_dtype)
        return rn.element(self.data.imag)

    @imag.setter
    def imag(self, newimag):
        """The setter for the imaginary part.

        This method is invoked by ``vec.imag = other``.

        Parameters
        ----------
        newreal : `array-like` or scalar
            The new imaginary part for this vector.

        Examples
        --------
        >>> x = Cn(3).element([5+1j, 3, 2-2j])
        >>> a = Rn(3).element([0, 0, 0])
        >>> x.imag = a; print(x)
        [(5+0j), (3+0j), (2+0j)]

        Other array-like types and broadcasting:

        >>> x.imag = 1.0; print(x)
        [(5+1j), (3+1j), (2+1j)]
        >>> x.imag = [0, 2, -1]; print(x)
        [(5+0j), (3+2j), (2-1j)]
        """
        self.imag.data[:] = newimag

    def conj(self, out=None):
        """The complex conjugate of this vector.

        Parameters
        ----------
        out : `FnVector`, optional
            Vector to which the complex conjugate is written.
            Must be an element of this vector's space.

        Returns
        -------
        out : `FnVector`
            The complex conjugate vector. If ``out`` was provided,
            the returned object is a reference to it.

        Examples
        --------
        >>> x = Cn(3).element([5+1j, 3, 2-2j])
        >>> y = x.conj(); print(y)
        [(5-1j), (3-0j), (2+2j)]

        The out parameter allows you to avoid a copy

        >>> z = Cn(3).element()
        >>> z_out = x.conj(out=z); print(z)
        [(5-1j), (3-0j), (2+2j)]
        >>> z_out is z
        True

        It can also be used for in-place conj

        >>> x_out = x.conj(out=x); print(x)
        [(5-1j), (3-0j), (2+2j)]
        >>> x_out is x
        True
        """
        if out is None:
            return self.space.element(self.data.conj())
        else:
            self.data.conj(out.data)
            return out

    def __ipow__(self, other):
        """Return ``self **= other``."""
        try:
            if other == int(other):
                return super().__ipow__(other)
        except TypeError:
            pass

        np.power(self.data, other, out=self.data)
        return self


if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
