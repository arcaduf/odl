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

"""Cone beam geometries in 3 dimensions."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

from odl.tomo.geometry.detector import Flat2dDetector
from odl.tomo.geometry.geometry import (
    DivergentBeamGeometry, AxisOrientedGeometry)

__all__ = ('FlexibleDivergentBeamGeometry',)


class FlexibleDivergentBeamGeometry(DivergentBeamGeometry):

    """An arbitrary divergent beam geometry

    See Also
    --------
    HelicalConeFlatGeometry
    CircularConeFlatGeometry
    """

    def __init__(self, apart, dpart, src_position, det_refpoint,
                 rotation_matrix, **kwargs):
        """Initialize a new instance.

        Parameters
        ----------
        apart : 1-dim. `RectPartition`
            Partition of the angle interval
        dpart : 2-dim. `RectPartition`
            Partition of the detector parameter rectangle


        Other Parameters
        ----------------

        """
        self._src_position = src_position
        self._det_refpoint = det_refpoint
        self._rotation_matrix = rotation_matrix

        assert callable(src_position)
        assert callable(det_refpoint)
        assert callable(rotation_matrix)

        det_init_axes = [[0, 1, 0], [0, 0, 1]]
        detector = Flat2dDetector(dpart, det_init_axes)

        super().__init__(ndim=3, motion_part=apart, detector=detector)

    @property
    def angles(self):
        """The discrete angles given in this geometry."""
        return self.motion_grid.coord_vectors[0]

    def det_refpoint(self, angle):
        """Return the detector reference point position at ``angle``.

        Parameters
        ----------
        angle : `float`
            Rotation angle given in radians, must be contained in
            this geometry's `motion_params`

        Returns
        -------
        point : `numpy.ndarray`, shape (3,)
            Detector reference point corresponding to the given angle

        See also
        --------
        rotation_matrix
        """
        angle = float(angle)
        if angle not in self.motion_params:
            raise ValueError('`angle` {} is not in the valid range {}'
                             ''.format(angle, self.motion_params))

        return self._det_refpoint(angle)

    def src_position(self, angle):
        """Return the source position at ``angle``.

        Parameters
        ----------
        angle : `float`
            Rotation angle given in radians, must be contained in
            this geometry's `motion_params`

        Returns
        -------
        point : `numpy.ndarray`, shape (3,)
            Detector reference point corresponding to the given angle

        See also
        --------
        rotation_matrix
        """
        angle = float(angle)
        if angle not in self.motion_params:
            raise ValueError('`angle` {} is not in the valid range {}'
                             ''.format(angle, self.motion_params))

        return self._src_position(angle)

    def rotation_matrix(self, angle):
        """The detector rotation function.

        Parameters
        ----------
        angle : `float`
            The motion parameter given in radian. It must be
            contained in this geometry's `motion_params`.

        Returns
        -------
        rot_mat : `numpy.ndarray`, shape ``(3, 3)``
            The rotation matrix mapping the standard basis vectors in
            the fixed ("lab") coordinate system to the basis vectors of
            the local coordinate system of the detector reference point,
            expressed in the fixed system.
        """
        angle = float(angle)
        if angle not in self.motion_params:
            raise ValueError('`angle` {} is not in the valid range {}'
                             ''.format(angle, self.motion_params))

        return self._rotation_matrix(angle)


if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
