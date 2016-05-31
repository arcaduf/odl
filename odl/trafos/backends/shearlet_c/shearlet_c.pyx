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

"""Cython wrapper for the Shearlet transform in 2 and 3 dimensions."""

import numpy as np
cimport numpy as np


__all__ = ('shearlet_trafo_2d',)

FLOAT = np.float32
ctypedef np.float32_t FLOAT_T
ctypedef FLOAT_T REAL  # TODO: use a macro to adjust this directive


cdef extern int shearlets2D(REAL *Image,
                            int Lx,
                            int Ly,
                            int wfilterlength,
                            REAL *wfilter,
                            int Level0,
                            int Level1,
                            int nsector,
                            REAL *coefficients)


def shearlet_trafo_2d(np.ndarray[FLOAT_T, ndim=2] image not None,
                      np.ndarray[FLOAT_T, ndim=2] coeffs_out not None,
                      np.ndarray[FLOAT_T, ndim=1] wfilter not None,
                      int lowpass_levels,
                      int shearing_levels,
                      int angles_per_cone):
    """Two-dimensional Shearlet transform wrapped by Cython.

    Signature::

        shearlet_trafo_2d(np.ndarray[float, ndim=2] image,
                      np.ndarray[float, ndim=2] coeffs_out,
                      np.ndarray[float, ndim=1] wfilter,
                      int lowpass_levels,
                      int shearing_levels,
                      int angles_per_cone)
    """
    assert (image.dtype == FLOAT and
            wfilter.dtype == FLOAT and
            coeffs_out.dtype == FLOAT)
    assert lowpass_levels > 0 and shearing_levels > 0 and angles_per_cone > 0

    shearlets2D(&image[0, 0], image.shape[0], image.shape[1],
                wfilter.shape[0], &wfilter[0],
                lowpass_levels, shearing_levels, angles_per_cone,
                &coeffs_out[0, 0])

