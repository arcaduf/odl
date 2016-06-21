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

"""Example using the ray transform with circular cone beam geometry."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

import numpy as np
import odl


# Discrete reconstruction space: discretized functions on the cube
# [-20, 20]^3 with 300 samples per dimension.
reco_space = odl.uniform_discr(
    min_corner=[-20, -20, -20], max_corner=[20, 20, 20],
    nsamples=[300, 300, 300], dtype='float32')

# Make a circular cone beam geometry with flat detector
# Angles: uniformly spaced, n = 360, min = 0, max = 2 * pi
angle_partition = odl.uniform_partition(-50, 50, 36)
# Detector: uniformly sampled, n = (558, 558), min = (-30, -30), max = (30, 30)
detector_partition = odl.uniform_partition([-50, -50], [50, 50], [558, 558])
#detector_partition = odl.uniform_partition([-100, -50], [100, 50], [558, 558])

d = 100
src_position = lambda x: np.array([-d + x, 0, 0])
det_refpoint = lambda x: np.array([50, 0, 0])
rotation_matrix = lambda x: np.eye(3)

geometry = odl.tomo.FlexibleDivergentBeamGeometry(
    angle_partition, detector_partition, src_position, det_refpoint,
    rotation_matrix)

# ray transform aka forward projection. We use ASTRA CUDA backend.
ray_trafo = odl.tomo.RayTransform(reco_space, geometry, impl='astra_cuda')

# Create a discrete Shepp-Logan phantom (modified version)
phantom = odl.util.phantom.shepp_logan(reco_space, True)

# Create projection data by calling the ray transform on the phantom
proj_data = ray_trafo(phantom)

# Back-projection can be done by simply calling the adjoint operator on the
# projection data (or any element in the projection space).
backproj = ray_trafo.adjoint(proj_data)

# Shows a slice of the phantom, projections, and reconstruction
phantom.show(coords=[None, 0, None], title='Phantom, middle y slice')
proj_data.show(coords=[0, None, None], title='Projection angle=0')
backproj.show(coords=[None, 0, None], title='back-projection, middle y slice')

print(proj_data.inner(proj_data))
print(phantom.inner(backproj))
print(proj_data.inner(proj_data) / phantom.inner(backproj))

#x = ray_trafo.domain.element()
#odl.solvers.conjugate_gradient_normal(ray_trafo, x, proj_data, 100, callback=odl.solvers.CallbackShow(coords=[None, 0, None]) & odl.solvers.CallbackPrintIteration())