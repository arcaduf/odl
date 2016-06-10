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
# along with ODL. If not, see <http://www.gnu.org/licenses/>.

# cython: profile=True

"""Cython implementation of the variable Lp proximal."""

import numpy as np

cimport cython
cimport numpy as np


__all__ = ('compute_varlp_prox_scalar',)


def compute_varlp_prox_scalar(f, p, out, s, max_newton_iter=5):
    """Compute the proximal of the variable Lp modular, scalar version.

    Signature::

        compute_varlp_prox_scalar(np.ndarray f,
                                  np.ndarray p,
                                  np.ndarray out,
                                  float s,
                                  int max_newton_iter)
    """
    if (f.dtype == np.float32 and
            p.dtype == np.float32 and
            out.dtype == np.float32):
        compute_varlp_prox_scalar_f32(f.ravel(), p.ravel(), out.ravel(),
                                      np.float32(s), max_newton_iter)
    elif (f.dtype == np.float64 and
          p.dtype == np.float64 and
          out.dtype == np.float64):
        compute_varlp_prox_scalar_f64(f.ravel(), p.ravel(), out.ravel(),
                                      np.float64(s), max_newton_iter)
    else:
        raise ValueError('no implementation for data types {}, {}, {} of '
                         '`f`, `p`, `out`'.format(f.dtype, p.dtype, out.dtype))


# --- Single precision --- #

FLOAT = np.float32
ctypedef np.float32_t FLOAT32_T


cdef extern void compute_varlp_prox_scalar__float(
    FLOAT32_T *f, FLOAT32_T *p, FLOAT32_T *out, int n, FLOAT32_T s,
    int max_newton_iter)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef compute_varlp_prox_scalar_f32(
        np.ndarray[FLOAT32_T] f,
        np.ndarray[FLOAT32_T] p,
        np.ndarray[FLOAT32_T] out,
        FLOAT32_T s,
        int max_newton_iter):
    """Compute the proximal of the variable Lp modular, scalar version.

    Signature::

        void compute_varlp_prox_scalar_f32(np.ndarray[float32] f,
                                           np.ndarray[float32] p,
                                           np.ndarray[float32] out,
                                           float32 s,
                                           int max_newton_iter)
    """
    assert (f.dtype == FLOAT and p.dtype == FLOAT and out.dtype == FLOAT)
    cdef int n = f.shape[0]
    compute_varlp_prox_scalar__float(&f[0], &p[0], &out[0], n, s,
                                     max_newton_iter)


# --- Double precision --- #

DOUBLE = np.float64
ctypedef np.float64_t FLOAT64_T


cdef extern void compute_varlp_prox_scalar__double(
    FLOAT64_T *f, FLOAT64_T *p, FLOAT64_T *out, int n, FLOAT64_T s,
    int max_newton_iter)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef compute_varlp_prox_scalar_f64(
        np.ndarray[FLOAT64_T] f,
        np.ndarray[FLOAT64_T] p,
        np.ndarray[FLOAT64_T] out,
        FLOAT64_T s,
        int max_newton_iter):
    """Compute the proximal of the variable Lp modular, scalar version.

    Signature::

        void compute_varlp_prox_scalar_f64(np.ndarray[float64] f,
                                           np.ndarray[float64] p,
                                           np.ndarray[float64] out,
                                           float64 s,
                                           int max_newton_iter)
    """
    assert (f.dtype == DOUBLE and p.dtype == DOUBLE and out.dtype == DOUBLE)
    cdef int n = f.shape[0]
    compute_varlp_prox_scalar__double(&f[0], &p[0], &out[0], n, s,
                                     max_newton_iter)
