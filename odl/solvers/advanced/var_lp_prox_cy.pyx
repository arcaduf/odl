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


__all__ = ('compute_varlp_prox_scalar_f32',)


FLOAT = np.float32
ctypedef np.float32_t FLOAT_T


cdef FLOAT_T startvalue(FLOAT_T val, FLOAT_T p, FLOAT_T sigma):
    cdef FLOAT_T sval = val / (p * (2.0 - p) * sigma)
    if sval <= 1:
        sval **= p - 1.0
    else:
        sval = 1.0
    return sval

@cython.profile(False)
cdef FLOAT_T numerator(FLOAT_T it, FLOAT_T p, FLOAT_T sigma, FLOAT_T val):
    return p * (p - 2) * sigma * it ** (p - 1) + val

@cython.profile(False)
cdef FLOAT_T denominator(FLOAT_T it, FLOAT_T p, FLOAT_T sigma):
    return 1.0 + p * (p - 1) * sigma * it ** (p - 2)

@cython.profile(False)
cdef FLOAT_T newton_iter_f32(FLOAT_T val, FLOAT_T sigma, FLOAT_T p,
                             int niter=5, FLOAT_T start_relax=0.5):
    """Helper function for the inner Newton iteration (single precision).

    Signature::

        float newton_iter_f32(float32 val, float32 sigma, float32 p,
                              int niter=5, float32 start_relax=0.5)
    """
    # Start value which guarantees convergence.
    cdef:
        int i
        FLOAT_T it, startval, denom, numer

#    startval = val / (p * (2.0 - p) * sigma)
#    if startval <= 1:
#        startval **= p - 1.0
#    else:
#        startval = 1.0
    startval = startvalue(val, p, sigma)

    startval *= start_relax

    # The iteration itself
    it = startval
    for i in range(niter):
        # Denominator 1 + p*(p-1)*sigma*q**(p-2)
#        denom = 1.0 + p * (p - 1) * sigma * it ** (p - 2)
        denom = denominator(it, p, sigma)

        # Numerator p*(p-2)*sigma*q**(p-1) + val
#        numer = p * (p - 2) * sigma * it ** (p - 1) + val
        numer = numerator(it, p, sigma, val)

        it = numer / denom

    return it


cdef inline FLOAT_T max(FLOAT_T x, FLOAT_T y):
    return x if x > y else y

cdef inline FLOAT_T min(FLOAT_T x, FLOAT_T y):
    return x if x < y else y

cdef inline FLOAT_T abs(FLOAT_T x):
    return x if x >= 0 else -x


cdef FLOAT_T varlp_prox_1(FLOAT_T f, FLOAT_T s):
    return max(1.0 - s / abs(f), 0.0) * f


cdef FLOAT_T varlp_prox_2(FLOAT_T f, FLOAT_T s):
    return f / (1.0 + 2.0 * s)


@cython.boundscheck(False)
@cython.wraparound(False)
def compute_varlp_prox_scalar_f32(
        np.ndarray[FLOAT_T] f,
        np.ndarray[FLOAT_T] p,
        np.ndarray[FLOAT_T] out,
        FLOAT_T s,
        int max_newton_iter):
    """Compute the proximal of the variable Lp modular, scalar version.

    Signature::

        void compute_varlp_prox_scalar_f32(np.ndarray[float32] f,
                                           np.ndarray[float32] p,
                                           np.ndarray[float32] out,
                                           float32 s,
                                           int max_newton_iter)
    """
    assert f.dtype == np.float32 and p.dtype == np.float32
    cdef:
        int i
        int nx = f.shape[0]

    for i in range(nx):
        if f[i] == 0.0:
            out[i] = 0.0

        elif p[i] <= 1.05:
#            out[i] = max(1.0 - s / abs(f[i]), 0.0) * f[i]
            out[i] = varlp_prox_1(f[i], s)

        elif p[i] >= 1.95:
#            out[i] = f[i] / (1.0 + 2.0 * s)
            out[i] = varlp_prox_2(f[i], s)

        else:
            out[i] = newton_iter_f32(f[i], s, p[i], max_newton_iter)

    return out

#
#def _call_scalar(self, f, out, **kwargs):
#    """Implement ``self(x, out, **kwargs)`` for scalar domain."""
#    if self.g is not None:
#        f = f - self.g
#
#    step = self.sigma * float(lam)
#
#    exp_arr = self.exponent.asarray()
#    out_arr = out.asarray()
#    f_arr = f.asarray()
#    f_nz = (f_arr != 0)
#
#    # p = 2
#    # This formula is used globally since it sets out to 0
#    # where f is 0.
#    out.lincomb(0, out, 1.0 / (1.0 + 2.0 * step), f)
#
#    # p = 1 (taking also close to one for stability)
#    cur_exp = (exp_arr <= 1.05)
#    current = cur_exp & f_nz
#    factor = np.maximum(1.0 - step / np.abs(f_arr[current]), 0.0)
#    out_arr[current] = factor * f_arr[current]
#
#    # Newton iteration for other p values. We consider only those
#    # entries that correspond to f != 0.
#    cur_exp = ~((exp_arr >= 1.95) | cur_exp)
#    current = cur_exp & f_nz
#    exp_p = exp_arr[current]
#    exp_m1 = exp_p - 1
#    exp_m2 = exp_p - 2
#    it = out_arr[current]
#    val = f_arr[current]
#    tmp = np.empty_like(it)
#
#    maxiter = int(kwargs.pop('max_newton_iter', 5))
#    self._newton_iter(it, np.abs(val), step, exp_p, exp_m1, exp_m2,
#                      niter=maxiter, tmp=tmp)
#
#    out_arr[current] = it
#    out_arr[current] *= np.sign(val)
#    out[:] = out_arr
#
#    if self.g is not None:
#        out += self.g
#
#def _call_pspace(self, f, out, **kwargs):
#    """Implement ``self(x, out, **kwargs)`` for vectorial domain."""
#    if self.g is not None:
#        f = f - self.g
#
#    step = self.sigma * float(lam)
#
#    exp_arr = self.exponent.asarray()
#    f_nz = [(fi.asarray() != 0) for fi in f]
#    pw_norm = PointwiseNorm(self.domain)
#    f_norm = pw_norm(f)
#
#    # p = 2
#    # This formula is used globally since it sets out to 0
#    # where f is 0.
#    out.lincomb(0, out, 1.0 / (1.0 + 2.0 * step), f)
#
#    # p = 1 (taking also close to one for stability)
#    cur_exp = (exp_arr <= 1.05)
#    for fi, fi_nz, oi in zip(f, f_nz, out):
#        fi_arr = fi.asarray()
#        oi_arr = oi.asarray()
#        current = cur_exp & fi_nz
#        factor = np.maximum(1.0 - step / np.abs(fi_arr[current]), 0.0)
#        oi_arr[current] = factor * fi_arr[current]
#        oi[:] = oi_arr
#
#    # Newton iteration for other p values
#    maxiter = int(kwargs.pop('max_newton_iter', 5))
#    cur_exp = ~((exp_arr >= 1.95) | cur_exp)
#    for fi, fi_nz, oi in zip(f, f_nz, out):
#        fi_arr = fi.asarray()
#        oi_arr = oi.asarray()
#        current = cur_exp & fi_nz
#        exp_p = exp_arr[current]
#        exp_m1 = exp_p - 1
#        exp_m2 = exp_p - 2
#
#        it = oi_arr[current]
#        val = fi_arr[current]
#        tmp = np.empty_like(it)
#        self._newton_iter(it, np.abs(val), step, exp_p, exp_m1, exp_m2,
#                          niter=maxiter, tmp=tmp)
#
#        oi_arr[current] = it
#        oi_arr[current] /= f_norm.asarray()[current]
#        oi_arr[current] *= val
#
#        oi[:] = oi_arr
#
#    if self.g is not None:
#        out += self.g
