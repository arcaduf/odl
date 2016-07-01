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

"""Example using the ray transform with 2d parallel beam geometry."""

import numpy as np
import odl
import time

# Discrete reconstruction space: discretized functions on the rectangle
# [-20, 20]^2 with 300 samples per dimension.
n  = 256
na = 2000
reco_space = odl.uniform_discr(
    min_corner=[-n,-n], max_corner=[n, n], nsamples=[n, n],
    dtype='float32')

# Make a parallel beam geometry with flat detector
# Angles: uniformly spaced, n = 360, min = 0, max = 2 * pi
angle_partition = odl.uniform_partition(0, np.pi, na)
# Detector: uniformly sampled, n = 558, min = -30, max = 30
detector_partition = odl.uniform_partition(-n, n, n)
geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)

# ray transform aka forward projection. We use 'scikit' backend.

#%%
ray_trafo = odl.tomo.RayTransform(reco_space, geometry, impl='fgp')


#%%
# Create a discrete Shepp-Logan phantom (modified version)
phantom = odl.util.phantom.shepp_logan(reco_space, modified=True)

#%%
# Create projection data by calling the ray transform on the phantom
proj_data = ray_trafo(phantom)
proj_data.show(title='Projection data (sinogram)')
#%%
# Back-projection can be done by simply calling the adjoint operator on the
# projection data (or any element in the projection space).
backproj = ray_trafo.adjoint(proj_data)

# Shows a slice of the phantom, projections, and reconstruction
phantom.show(title='Phantom')
proj_data.show(title='Projection data (sinogram)')
backproj.show(title='Back-projected data')


'''
# Test speed for the computation of the forward projection
n  = 128 #2048
na = np.array( [ 800 , 1600 , 2400 , 3200 ] )

for i in range( len( na ) ):
    reco_space = odl.uniform_discr(
        min_corner=[-20, -20], max_corner=[20, 20], nsamples=[n, n],
        dtype='float32')    
    angle_partition = odl.uniform_partition(0, np.pi, na[i])
    detector_partition = odl.uniform_partition(-30, 30, n)
    geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)
    phantom = odl.util.phantom.shepp_logan(reco_space, modified=True)
    ray_trafo = odl.tomo.RayTransform(reco_space, geometry, impl='fgp')
    time1 = time.time()
    proj_data = ray_trafo(phantom)
    time2 = time.time()
    print( '\nTime elapsed to compute sinogram ', n,' pixels X ', na[i], ' views: ' , time2 - time1 )
'''
