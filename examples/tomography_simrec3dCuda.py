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

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

from math import sin, cos, pi
import matplotlib.pyplot as plt
import numpy as np

import odl
import SimRec2DPy as SR


class ProjectionGeometry3D(object):
    """ Geometry for a specific projection
    """
    def __init__(self, sourcePosition, detectorOrigin, pixelDirectionU,
                 pixelDirectionV):
        self.sourcePosition = sourcePosition
        self.detectorOrigin = detectorOrigin
        self.pixelDirectionU = pixelDirectionU
        self.pixelDirectionV = pixelDirectionV

    def __repr__(self):
        return ('ProjectionGeometry3D({}, {}, {}, {})'
                ''.format(self.sourcePosition, self.detectorOrigin,
                          self.pixelDirectionU, self.pixelDirectionV))


class CudaProjector3D(odl.LinearOperator):
    """ A projector that creates several projections as defined by geometries
    """
    def __init__(self, volumeOrigin, voxelSize, nVoxels, nPixels, stepSize,
                 geometries, domain, range):
        super().__init__(domain, range)
        self.geometries = geometries
        self.forward = SR.SRPyCuda.CudaForwardProjector3D(
            nVoxels, volumeOrigin, voxelSize, nPixels, stepSize)
        self._adjoint = CudaBackProjector3D(
            volumeOrigin, voxelSize, nVoxels, nPixels, stepSize, geometries,
            range, domain)

    @odl.util.timeit("Project")
    def _apply(self, volume, projection):
        # Create projector
        self.forward.setData(volume.ntuple.data_ptr)

        # Project all geometries
        for i in range(len(self.geometries)):
            geo = self.geometries[i]

            self.forward.project(geo.sourcePosition, geo.detectorOrigin,
                                 geo.pixelDirectionU, geo.pixelDirectionV,
                                 projection[i].ntuple.data_ptr)

    @property
    def adjoint(self):
        return self._adjoint


class CudaBackProjector3D(odl.LinearOperator):
    def __init__(self, volumeOrigin, voxelSize, nVoxels, nPixels, stepSize,
                 geometries, domain, range):
        super().__init__(domain, range)
        self.geometries = geometries
        self.back = SR.SRPyCuda.CudaBackProjector3D(nVoxels, volumeOrigin,
                                                    voxelSize, nPixels,
                                                    stepSize)

    @odl.util.timeit("BackProject")
    def _apply(self, projections, out):
        # Zero out the return data
        out.set_zero()

        # Append all projections
        for geo, proj in zip(self.geometries, projections):
            self.back.backProject(
                geo.sourcePosition, geo.detectorOrigin, geo.pixelDirectionU,
                geo.pixelDirectionV, proj.ntuple.data_ptr, out.ntuple.data_ptr)

        # correct for unmatched projectors
        out *= 3.0  # 0.165160099353932

# Set geometry parameters
volumeSize = np.array([224.0, 224.0, 1.0])
volumeOrigin = np.array([-112.0, -112.0, 10.0])  # -volumeSize/2.0

detectorSize = np.array([287.04, 264.94])
detectorOrigin = np.array([-143.52, 0.0])

sourceAxisDistance = 790.0
detectorAxisDistance = 210.0

# Discretization parameters
# nVoxels, nPixels = np.array([44, 44, 27]), np.array([78, 72])
nVoxels, nPixels = np.array([448, 448, 1]), np.array([780, 720])
nProjection = 332

# Scale factors
voxelSize = volumeSize/nVoxels
pixelSize = detectorSize/nPixels
stepSize = voxelSize.max()/1.0

# Define projection geometries
geometries = []
for theta in np.linspace(0, 2*pi, nProjection, endpoint=False):
    x0 = np.array([cos(theta), sin(theta), 0.0])
    y0 = np.array([-sin(theta), cos(theta), 0.0])
    z0 = np.array([0.0, 0.0, 1.0])

    projSourcePosition = -sourceAxisDistance * x0
    projPixelDirectionU = y0 * pixelSize[0]
    projPixelDirectionV = z0 * pixelSize[1]
    projDetectorOrigin = (detectorAxisDistance * x0 + detectorOrigin[0] * y0 +
                          detectorOrigin[1] * z0)
    geometries.append(ProjectionGeometry3D(
        projSourcePosition, projDetectorOrigin, projPixelDirectionU,
        projPixelDirectionV))

# Define the space of one projection
projectionSpace = odl.L2(odl.Rectangle([0, 0], detectorSize))

# Discretize projection space
projectionDisc = odl.l2_uniform_discretization(projectionSpace, nPixels,
                                               impl='cuda', order='F')

# Create the data space, which is the Cartesian product of the
# single projection spaces
dataDisc = odl.ProductSpace(projectionDisc, nProjection)

# Define the reconstruction space
reconSpace = odl.L2(odl.Cuboid([0, 0, 0], volumeSize))

# Discretize the reconstruction space
reconDisc = odl.l2_uniform_discretization(reconSpace, nVoxels,
                                          impl='cuda', order='F')

# Create a phantom
phantom = SR.SRPyUtils.phantom(nVoxels[0:2])
phantom = np.repeat(phantom, nVoxels[-1]).reshape(nVoxels)
phantomVec = reconDisc.element(phantom)

# Make the operator
projector = CudaProjector3D(volumeOrigin, voxelSize, nVoxels, nPixels,
                            stepSize, geometries, reconDisc, dataDisc)

projections = projector(phantomVec)
"""
#Apply once to find norm estimate
recon = projector.T(projections)
#print('normEst', recon.norm() / phantomVec.norm())
print('const = ', projections.inner(projector(phantomVec))/ phantomVec.inner(projector.adjoint(projections)))
raise Exception()
del recon"""
del phantomVec


# Define function to plot each result
def plotResult(x):
    with odl.util.Timer('plotting'):
        plt.figure()
        plt.imshow(x.asarray()[:, :, nVoxels[2]//2])
        #plt.clim(0, 1)
        plt.colorbar()
        plt.show()

# Solve using landweber
x = reconDisc.zero()
#odl.operator.solvers.landweber(projector, x, projections, 100, omega=0.5/normEst,
#                               partial=odl.operator.solvers.ForEachPartial(plotResult))
# solvers.landweber(projector, x, projections, 10, omega=0.4/normEst,
#                   partial=solvers.PrintIterationPartial())
odl.operator.solvers.conjugate_gradient_normal(
    projector, x, projections, 100,
    partial=odl.operator.solvers.ForEachPartial(plotResult))