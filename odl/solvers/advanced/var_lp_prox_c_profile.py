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

import cProfile
import numpy as np
import pstats

import odl

import pyximport
pyximport.install()

from var_lp_prox_c import compute_varlp_prox_scalar_f32

space = odl.uniform_discr([-1, -1], [1, 1], (1024, 1024), dtype='float32')

f = 2.0 * space.one()
f_arr = f.ntuple.data

p_randarr = np.random.uniform(low=1, high=2, size=space.shape)
p = space.element(p_randarr)

out_arr = np.empty_like(f_arr)

prox_op = odl.solvers.advanced.proximal_variable_lp(space, p)(0.1)
p_arr = p.ntuple.data

cProfile.runctx('compute_varlp_prox_scalar_f32(f_arr, p_arr, out_arr, 0.1, 3)',
                globals(), locals(), 'varlp_prox.prof')

s = pstats.Stats('varlp_prox.prof')
s.strip_dirs().sort_stats('time').print_stats()
