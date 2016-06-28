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

"""CPU implementations of ``n``-dimensional Cartesian spaces."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

import numpy as np

from odl.set.sets import Set, RealNumbers, ComplexNumbers
from odl.set.space import LinearSpace, LinearSpaceVector
from odl.space.base_ntuples import _TYPE_MAP_R2C, _TYPE_MAP_C2R
from odl.util.ufuncs import NtuplesUFuncs
from odl.util.utility import (
    arraynd_repr, dtype_repr, is_scalar_dtype, is_floating_dtype,
    is_real_dtype, is_complex_floating_dtype)


__all__ = ('TensorSetBase', 'PreBaseTensor', 'TensorSpaceBase', 'BaseTensor')


class TensorSetBase(Set):

    """Base class for sets of tensors of arbitrary type."""

    def __init__(self, shape, dtype, order):
        """Initialize a new instance.

        Parameters
        ----------
        shape : sequence of int
            Number of elements per axis.
        dtype :
            Data type of each element. Can be provided in any
            way the `numpy.dtype` function understands, most notably
            as built-in type, as one of NumPy's internal datatype
            objects or as string.
        order : {'C', 'F'}
            Axis ordering of the data storage.
        """
        self._shape = tuple(int(s) for s in shape)
        if any(s < 0 for s in self.shape):
            raise ValueError('`shape` must have only positive entries, got '
                             '{}'.format(shape))
        self._dtype = np.dtype(dtype)
        self._order = str(order).upper()
        if self.order not in ('C', 'F'):
            raise ValueError("`order '' not understood".format(order))

    @property
    def shape(self):
        """Number of elements per axis."""
        return self._shape

    @property
    def dtype(self):
        """Data type of each tensor entry."""
        return self._dtype

    @property
    def order(self):
        """Data storage order, either 'C' or 'F'."""
        return self._order

    @property
    def size(self):
        """Total number of tensor entries."""
        return np.prod(self.shape)

    @property
    def ndim(self):
        """Number of dimensions (axes)."""
        return len(self.shape)

    def __len__(self):
        """Number of tensor entries along the first axis."""
        return self.shape[0]

    def __contains__(self, other):
        """Return ``other in self``.

        Returns
        -------
        contains : bool
            True if ``other.space`` is equal to ``self``,
            False otherwise.

        Examples
        --------
        TODO: fix example

        >>> import odl
        >>> mat_space = NumpyTensors([2, 3], dtype='float32')
        >>> long_3 = Ntuples(3, dtype='int64')
        >>> long_3.element() in long_3

        >>> long_3.element() in Ntuples(3, dtype='int32')

        >>> long_3.element() in Ntuples(3, dtype='float64')

        """
        return getattr(other, 'space', None) == self

    def __eq__(self, other):
        """Return ``self == other``.

        Returns
        -------
        equals : bool
            True if ``other`` is an instance of ``type(self)`` with the
            same `shape`, `dtype` and `order`, otherwise `False`.

        Examples
        --------
        TODO: example

        >>> from odl import Ntuples
        >>> int_3 = Ntuples(3, dtype=int)
        >>> int_3 == int_3


        Equality is not identity:

        >>> int_3a, int_3b = Ntuples(3, int), Ntuples(3, int)
        >>> int_3a == int_3b

        >>> int_3a is int_3b


        >>> int_3, int_4 = Ntuples(3, int), Ntuples(4, int)
        >>> int_3 == int_4

        >>> int_3, str_3 = Ntuples(3, 'int'), Ntuples(3, 'S2')
        >>> int_3 == str_3

        """
        if other is self:
            return True

        return ((isinstance(self, type(other)) or
                 isinstance(other, type(self))) and
                self.shape == other.shape and
                self.dtype == other.dtype and
                self.order == other.order)

    def __repr__(self):
        """Return ``repr(self)``."""
        return "{}({}, {}, '{}')".format(self.__class__.__name__, self.shape,
                                         dtype_repr(self.dtype), self.order)

    @property
    def element_type(self):
        """Type of elements in this set: `PreBaseTensor`."""
        return PreBaseTensor


class PreBaseTensor(object):

    """Abstract class for representation of `TensorSetBase` elements.

    Defines abstract and concrete attributes independent of data
    representation.
    """

    def __init__(self, space, *args, **kwargs):
        """Initialize a new instance."""
        self._space = space

    def copy(self):
        """Create an identical (deep) copy of this tensor."""
        raise NotImplementedError

    def asarray(self, out=None):
        """Extract the data of this tensor as a Numpy array.

        Parameters
        ----------
        out : `numpy.ndarray`
            Array to write the result to.

        Returns
        -------
        asarray : `numpy.ndarray`
            Numpy array of the same data type and shape as the space.
            If ``out`` was given, the returned object is a reference
            to it.
        """
        raise NotImplementedError

    def __getitem__(self, indices):
        """Return ``self[indices]``.

        Parameters
        ----------
        indices : index expression
            Integer, slice or sequence of these, defining the positions
            of the data array which should be accessed.

        Returns
        -------
        values : `TensorSetBase.dtype` or `PreBaseTensor`
            The value(s) at the given indices. Note that depending on
            the implementation, the returned object may be a (writable)
            view into the original array.
        """
        raise NotImplementedError

    def __setitem__(self, indices, values):
        """Implement ``self[indices] = values``.

        Parameters
        ----------
        indices : index expression
            Integer, slice or sequence of these, defining the positions
            of the data array which should be written to.
        values : scalar, array-like or `PreBaseTensor`
            The value(s) that are to be assigned.

            If ``index`` is an integer, ``value`` must be a scalar.

            If ``index`` is a slice or a sequence of slices, ``value``
            must be broadcastable to the shape of the slice.
        """
        raise NotImplementedError

    def __eq__(self, other):
        """Return ``self == other``.

        Returns
        -------
        equals : bool
            True if all entries of ``self`` and ``other`` are equal,
            False otherwise.
        """
        raise NotImplementedError

    @property
    def space(self):
        """Tensor space in which this tensor lives."""
        return self._space

    @property
    def shape(self):
        """Number of elements per axis."""
        return self.space.shape

    @property
    def dtype(self):
        """Data type of each entry."""
        return self.space.dtype

    @property
    def order(self):
        """Data storage order, either 'C' or 'F'."""
        return self._order

    @property
    def size(self):
        """Total number of entries."""
        return self.space.size

    @property
    def ndim(self):
        """Number of dimensions (axes)."""
        return self.space.ndim

    def __len__(self):
        """Return ``len(self)``.

        Return the number of space dimensions.
        """
        return len(self.space)

    @property
    def itemsize(self):
        """Size in bytes of one tensor entry."""
        return self.dtype.itemsize

    @property
    def nbytes(self):
        """Number of bytes in memory occupied by this tensor."""
        return self.size * self.itemsize

    def __array__(self, dtype=None):
        """Return a Numpy array from this tensor.

        Parameters
        ----------
        dtype :
            Specifier for the data type of the output array.

        Returns
        -------
        array : `numpy.ndarray`
        """
        if dtype is None:
            return self.asarray()
        else:
            return self.asarray().astype(dtype, copy=False)

    def __array_wrap__(self, obj):
        """Return a new tensor wrapping the array ``obj``.

        Parameters
        ----------
        obj : `numpy.ndarray`
            Array to be wrapped.

        Returns
        -------
        vector : `PreBaseTensor`
            Tensor wrapping ``obj``.
        """
        if obj.ndim == 0:
            return self.space.field.element(obj)
        else:
            return self.space.element(obj)

    def __ne__(self, other):
        """Return ``self != other``."""
        return not self.__eq__(other)

    def __str__(self):
        """Return ``str(self)``."""
        return arraynd_repr(self)

    def __repr__(self):
        """Return ``repr(self)``."""
        return '{!r}.element({})'.format(self.space,
                                         arraynd_repr(self))

    @property
    def ufunc(self):
        """`TensorSetBaseUfunc`, access to numpy style ufuncs.

        These are always available, but may or may not be optimized for
        the specific space in use.
        """
        # TODO: add!
        raise NotImplementedError

    def show(self, title=None, method='scatter', show=False, fig=None,
             **kwargs):
        """Display this tensor graphically for ``ndim == 1 or 2``.

        Parameters
        ----------
        title : `str`, optional
            Set the title of the figure

        method : `str`, optional
            1d methods:

            'plot' : graph plot

            'scatter' : point plot

        show : `bool`, optional
            If the plot should be showed now or deferred until later.

        fig : `matplotlib.figure.Figure`
            The figure to show in. Expected to be of same "style", as
            the figure given by this function. The most common use case
            is that ``fig`` is the return value from an earlier call to
            this function.

        kwargs : {'figsize', 'saveto', ...}
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
        # TODO: take code from DiscreteLpVector and adapt
        from odl.util.graphics import show_discrete_data
        from odl.discr import RegularGrid
        grid = RegularGrid(0, self.size - 1, self.size)
        return show_discrete_data(self.asarray(), grid, title=title,
                                  method=method, show=show, fig=fig, **kwargs)


class TensorSpaceBase(TensorSetBase, LinearSpace):

    """Base class for tensor spaces independent of implementation."""

    def __init__(self, shape, dtype, order):
        """Initialize a new instance.

        Parameters
        ----------
        shape : sequence of int
            Number of elements per axis.
        dtype :
            Data type of each element. Can be provided in any
            way the `numpy.dtype` function understands, most notably
            as built-in type, as one of NumPy's internal datatype
            objects or as string.
            Only scalar data types (numbers) are allowed.
        order : {'C', 'F'}
            Axis ordering of the data storage.
        """
        TensorSetBase.__init__(self, shape, dtype, order)

        if not is_scalar_dtype(self.dtype):
            raise TypeError('`dtype` must be a scalar data type, got {!r}'
                            ''.format(dtype))

        if is_real_dtype(self.dtype):
            field = RealNumbers()
            self._is_real = True
            self._real_dtype = self.dtype
            self._real_space = self
            self._complex_dtype = _TYPE_MAP_R2C.get(self.dtype, None)
            self._complex_space = None  # Set in first call of astype
        else:
            field = ComplexNumbers()
            self._is_real = False
            self._real_dtype = _TYPE_MAP_C2R[self.dtype]
            self._real_space = None  # Set in first call of astype
            self._complex_dtype = self.dtype
            self._complex_space = self

        self._is_floating = is_floating_dtype(self.dtype)
        LinearSpace.__init__(self, field)

    @property
    def is_real_tensor_space(self):
        """True if this is a space of real tensors."""
        return self._is_real and self._is_floating

    @property
    def is_complex_tensor_space(self):
        """True if this is a space of complex tensors."""
        return (not self._is_real) and self._is_floating

    def _astype(self, dtype):
        """Internal helper for ``astype``."""
        # TODO: weighting
        return type(self)(self.shape, dtype=dtype, order=self.order)

    def astype(self, dtype):
        """Return a copy of this space with new ``dtype``.

        Parameters
        ----------
        dtype :
            Data type of the returned space. Can be given in any way
            `numpy.dtype` understands, e.g. as string (``'complex64'``)
            or data type (``complex``).

        Returns
        -------
        newspace : `TensorSpaceBase`
            Version of this space with given data type.
        """
        if dtype is None:
            # Need to filter this out since Numpy iterprets it as 'float'
            raise ValueError('unknown data type `None`')

        dtype = np.dtype(dtype)
        if dtype == self.dtype:
            return self

        # Caching for real and complex versions (exact dtyoe mappings)
        if dtype == self._real_dtype:
            if self._real_space is None:
                self._real_space = self._astype(dtype)
            return self._real_space
        elif dtype == self._complex_dtype:
            if self._complex_space is None:
                self._complex_space = self._astype(dtype)
            return self._complex_space
        else:
            return self._astype(dtype)

    @property
    def examples(self):
        """Return example random vectors."""
        # Always return the same numbers
        rand_state = np.random.get_state()
        np.random.seed(1337)

        yield ('Linspaced', self.element(
            np.linspace(0, 1, self.size).reshape(self.shape)))

        if self.is_real_tensor_space:
            yield ('Random noise', self.element(np.random.rand(*self.shape)))
        elif self.is_complex_tensor_space:
            yield ('Random noise',
                   self.element(np.random.rand(*self.shape) +
                                np.random.rand(*self.shape) * 1j))

        yield ('Normally distributed random noise',
               self.element(np.random.randn(self.size)))

        np.random.set_state(rand_state)

    def zero(self):
        """Return a tensor of all zeros."""
        raise NotImplementedError

    def one(self):
        """Return a tensor of all ones."""
        raise NotImplementedError

    def _multiply(self, x1, x2, out):
        """The entry-wise product of two tensors, assigned to ``out``."""
        raise NotImplementedError

    def _divide(self, x1, x2, out):
        """The entry-wise quotient of two tensors, assigned to ``out``."""
        raise NotImplementedError

    @property
    def element_type(self):
        """Type of elements in this set: `BaseTensor`."""
        return BaseTensor


class BaseTensor(PreBaseTensor, LinearSpaceVector):

    """Abstract class for representation of `TensorSpaceBase` elements.

    Defines abstract and concrete attributes and methods independent
    of data representation.
    """

    def __eq__(self, other):
        return LinearSpaceVector.__eq__(self, other)

    def copy(self):
        return LinearSpaceVector.copy(self)


if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
