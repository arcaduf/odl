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

"""Backends for other libraries."""

from __future__ import absolute_import

__all__ = ()

__all__ += ('stir_bindings',)

from . import stir_bindings

try:
    from . import astra_cpu
    from .astra_cpu import *
    __all__ += astra_cpu.__all__
    ASTRA_AVAILABLE = True
except ImportError:
    ASTRA_AVAILABLE = False

try:
    from . import astra_cuda
    from .astra_cuda import *
    __all__ += astra_cuda.__all__
    ASTRA_CUDA_AVAILABLE = True
except ImportError:
    ASTRA_CUDA_AVAILABLE = False

if ASTRA_AVAILABLE:
    from . import astra_setup
    from .astra_setup import *
    __all__ += astra_setup.__all__

__all__ += ('ASTRA_AVAILABLE', 'ASTRA_CUDA_AVAILABLE')
