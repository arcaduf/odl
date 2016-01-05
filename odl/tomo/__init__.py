# Copyright 2014, 2015 The ODL development group
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

"""Tomography related operators and geometries based on ODL."""


from __future__ import absolute_import

__version__ = '0.9b1'
__all__ = ('backends', 'geometry', 'operators')


from . import backends
from .backends import *
__all__ += backends.__all__

from . import geometry
from .geometry import *
__all__ += geometry.__all__

from . import operators
from .operators import *
__all__ += operators.__all__
