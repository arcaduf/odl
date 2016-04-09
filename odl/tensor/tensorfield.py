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

"""Vector and Tensor fields."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from builtins import super, zip
from future import standard_library
from future.utils import raise_from
standard_library.install_aliases()

import numpy as np

from odl.discr.discretization import Discretization
from odl.space.fspace import FunctionSpace
from odl.space.pspace import ProductSpace


__all__ = ()


class DiscreteTensorField(ProductSpace):

    """Implementation of discretized tensor fields on ``R^d``.

    Tensor fields are power spaces of discretized function spaces
    with a given shape. For example, a ``(3, 3)`` tensor field can be
    considered a ``(3, 3)`` array of scalar function spaces. An element
    of the tensor space is in each component an element of the underlying
    scalar space.
    """

    def __init__(self, space, shape):
        """Initialize a new instance.

        Parameters
        ----------
        space : `Discretization`
            Base scalar space for each component of the tensor field
        shape : positive `int` or `sequence` of positive `int`
            Components of the tensor field space. For vector fields
            (rank-1 tensor fields), a single integer can be used.
        """
        if not isinstance(space, Discretization):
            raise TypeError('space {!r} is not a Discretization instance.'
                            ''.format(space))

        if not isinstance(space.uspace, FunctionSpace):
            raise TypeError('non-discretized space is not a FunctionSpace '
                            'instance.')

        size = np.size(shape)
        super().__init__(space, size)

        self._scalar_space = space

        if np.isscalar(shape):
            self._shape = (int(shape),)
        else:
            self._shape = tuple(int(s) for s in shape)

        if not np.all(np.greater(self.shape, 0)):
            raise ValueError('shape {} is not positive.'.format(shape))

        self._dtype = np.dtype((self.scalar_space.dtype, self.shape))

    @property
    def rank(self):
        """Tensor field rank, i.e. the number of tensor axes."""
        return len(self)

    @property
    def dtype(self):
        """Data type of the tensor field, including the shape.

        The data type uses Numpy's advanced data type functionality
        as ``dtype = np.dtype((scalar_dtype, shape))``.
        """
        return self._dtype

    @property
    def scalar_space(self):
        """The scalar base space of this tensor field space."""
        return self._scalar_space

if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
