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

"""
Example of shape-based image reconstruction
using optimal information transformation.
"""

# Initial setup
import numpy as np
import numba
import matplotlib.pyplot as plt
import time
import ddmatch
from odl.operator.operator import Operator
import odl


class DeformationOperator(Operator):
    """Operator mapping parameter to a fixed deformed template.

    This operator computes for a fixed template ``I`` the deformed
    template::

        invphi(.) --> I(invphi(.))

    where ``invphi`` is the deformation parameter as follows::

        invphi: x --> invphi(x)

    Here, ``x`` is an element in the domain of target (ground truth).
    """

    def __init__(self, template):
        """Initialize a new instance.

        Parameters
        ----------
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        """
        # The operator maps from the parameter space to the template
        # (image) space.
        self.template = template

        # Create the space for the general inverse deformation
        self.domain_space = odl.ProductSpace(self.template.space,
                                             self.template.space.ndim)

        super().__init__(self.domain_space, self.template.space, linear=False)

    def _call(self, invdeformation):
        """Implementation of ``self(invdeformation)``.

        Parameters
        ----------
        invdeformation: `ProductSpaceVector`
            General inverse deformation for image grid points.
        """
        for i in xrange(s):
            for j in xrange(s):
                xind = int(xphi[i, j])
                yind = int(yphi[i, j])
                xindp1 = xind+1
                yindp1 = yind+1
                deltax = xphi[i, j]-float(xind)
                deltay = yphi[i, j]-float(yind)

                # Id xdelta is negative it means that xphi is negative, so xind
                # is larger than xphi. We then reduce xind and xindp1 by 1 and
                # after that impose the periodic boundary conditions.
                if (deltax < 0 or xind < 0):
                    deltax += 1.0
                    xind -= 1
                    xind %= s
                    xindp1 -= 1
                    xindp1 %= s
                elif (xind >= s):
                    xind %= s
                    xindp1 %= s
                elif (xindp1 >= s):
                    xindp1 %= s

                if (deltay < 0 or xind < 0):
                    deltay += 1.0
                    yind -= 1
                    yind %= s
                    yindp1 -= 1
                    yindp1 %= s
                elif (yind >= s):
                    yind %= s
                    yindp1 %= s
                elif (yindp1 >= s):
                    yindp1 %= s

                onemdeltax = 1.-deltax
                onemdeltay = 1.-deltay
                Iout[i, j] = I[yind, xind]*onemdeltax*onemdeltay + \
                    I[yind, xindp1]*deltax*onemdeltay + \
                    I[yindp1, xind]*deltay*onemdeltax + \
                    I[yindp1, xindp1]*deltay*deltax

        image_pts = self.template.space.grid.points()
        image_pts += np.asarray(displacement).T

        return self.template.interpolation(image_pts.T, bounds_check=False)

    def linear_deform(self, template, displacement):
        """Implementation of ``self(template, displacement)``.

        Parameters
        ----------
        template : `DiscreteLpVector`
            Fixed template deformed by the vector field. Its space
            must have the same number of dimensions as ``par_space``.
        displacement: `ProductSpaceVector`
            Linearized deformation parameters for image grid points.
        """
        image_pts = template.space.grid.points()
        image_pts += np.asarray(displacement).T
        return template.interpolation(image_pts.T, bounds_check=False)

    def derivative(self, displacement):
        """Frechet derivative of this operator in ``disp``.

        Parameters
        ----------
        displacement: `ProductSpaceVector`
            Linearized deformation parameters for image grid points.

        Returns
        -------
        deriv : `Operator`
        The derivative of this operator, evaluated at ``displacement``
        """
        deriv_op = LinearizedDeformationDerivative(self.template,
                                                   displacement)
        return deriv_op


def generate_optimized_density_match_L2_gradient_rec(image):
    s = image.shape[0]
    if (len(image.shape) != 2):
        raise(NotImplementedError('Only 2d images are allowed so far.'))
    if (image.shape[1] != s):
        raise(NotImplementedError('Only square images are allowed so far.'))
    if (image.dtype != np.float64):
        raise(NotImplementedError('Only float64 images are allowed so far.'))

    @numba.njit('void(f8,f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:])')
    def density_match_L2_gradient_2d_rec(sigma, dsqrtJdx, dsqrtJdy,
                                         dtmpdx, dtmpdy, doutdx, doutdy):
        for i in xrange(s):
            for j in xrange(s):
                doutdx[i, j] = sigma*dsqrtJdx[i, j] + 2. * dtmpdx[i, j]
                doutdy[i, j] = sigma*dsqrtJdy[i, j] + 2. * dtmpdy[i, j]

    return density_match_L2_gradient_2d_rec


def proj_noise(proj_data_shape, mu=0.0, sigma=0.1):
    """Produce white Gaussian noise for projections of n-D images.

       Produce white Gaussian noise for projections, with the same size
       as the number of projections.

       Parameters
       ----------
       proj_data_shape : shape of the projections
           Give the size of noise
       mu : Mean of the norm distribution
           The defalt is 0.0
       sigma : Standard deviation of the norm distribution
           The defalt is 0.1.
    """

    return np.random.normal(mu, sigma, proj_data_shape)


def SNR(signal, noise):
    """Compute the signal-to-noise ratio in dB.
    This compute::

    SNR = 10 * log10 (
        |signal - mean(signal)| / |noise - mean(noise)|)

    Parameters
    ----------
    signal : projection
    noise : white noise
    """
    if np.abs(np.asarray(noise)).sum() != 0:
        ave1 = np.sum(signal)/signal.size
        ave2 = np.sum(noise)/noise.size
        en1 = np.sqrt(np.sum((signal - ave1) * (signal - ave1)))
        en2 = np.sqrt(np.sum((noise - ave2) * (noise - ave2)))

        return 10.0 * np.log10(en1/en2)
    else:
        return 10000.


I0name = '../ddmatch/Example3 letters/c_highres.png'
I1name = '../ddmatch/Example3 letters/i_highres.png'
# I0name = 'Example3 letters/eight.png'
# I1name = 'Example3 letters/b.png'
# I0name = 'Example3 letters/v.png'
# I1name = 'Example3 letters/j.png'
# I0name = 'Example9 letters big/V.png'
# I1name = 'Example9 letters big/J.png'
# I0name = 'Example11 skulls/handnew1.png'
# I1name = 'Example11 skulls/handnew2.png'
# I0name = 'Example8 brains/DS0002AxialSlice80.png'
# I1name = 'Example8 brains/DS0003AxialSlice80.png'

I0 = plt.imread(I0name).astype('float')
I1 = plt.imread(I1name).astype('float')

I0 = I0[::2, ::2]
I1 = I1[::2, ::2]

# Create 2-D discretization reconstruction space
# The size of the domain should be proportional to the given images
discr_space = odl.uniform_discr([-16, -16],
                                [16, 16], [128, 128],
                                dtype='float32', interp='linear')

# Create the ground truth as the given image
ground_truth = discr_space.element(I0.T)

# Create the template as the given image
template = discr_space.element(I1.T)

# Compose mass-preserving operator to template
template_mass_pre = template

# Give the number of directions
num_angles = 6

# Create the uniformly distributed directions
angle_partition = odl.uniform_partition(
    0, np.pi, num_angles, nodes_on_bdry=[(True, False)])

# Create 2-D projection domain
detector_partition = odl.uniform_partition(-24, 24, 192)

# Create 2-D parallel projection geometry
geometry = odl.tomo.Parallel2dGeometry(angle_partition,
                                       detector_partition)

# Create projection data by given setting
xray_trafo_op = odl.tomo.RayTransform(discr_space, geometry, impl='astra_cuda')

# Create projection data by given setting
proj_data = xray_trafo_op(ground_truth)

# Create white Gaussian noise
noise = 5.0 * proj_data.space.element(proj_noise(proj_data.shape))

# Output the signal-to-noise ratio
print('snr = {!r}'.format(SNR(proj_data, noise)))

# Create noisy projection data
noise_proj_data = proj_data + noise

# Do the backprojection reconstruction
backproj = xray_trafo_op.adjoint(noise_proj_data)

density_match_L2_gradient_rec = \
    generate_optimized_density_match_L2_gradient_rec(I1)

# Regularization parameter, should be nonnegtive
sigma = 1000e-1

# Step size for the gradient descent method
epsilon = 0.005

# Maximum iteration number
n_iter = 500

dm = ddmatch.TwoComponentDensityMatching(source=I1, target=I0, sigma=sigma)

# Normalize the mass of template as ground truth
W1 = I1.T * np.linalg.norm(I0, 'fro')/np.linalg.norm(I1, 'fro')

# Normalized template as an element of discretization space
template = discr_space.element(W1)

# Create the memory for energy in each iteration
E = []
kE = len(E)
E = np.hstack((E, np.zeros(n_iter)))

axis = [np.linspace(
    0, I1.shape[i], I1.shape[i], endpoint=False) for i in range(I1.ndim)]
id_map = np.meshgrid(*axis)

vec_field = [np.zeros_like(I1) for _ in range(I1.ndim)]

vx = np.zeros_like(I1)
vy = np.zeros_like(I1)

tmp = list(id_map)

tmpx = dm.idx.copy()
tmpy = dm.idy.copy()

# Test time, set starting time
start = time.clock()

for k in xrange(n_iter):

    # OUTPUT
    E[k+kE] = (sigma*(dm.sqrtJ - 1)**2).sum()

    # STEP 1: update template_mass_pre
    template_array = np.asarray(template, dtype='float64')
    template_mass_pre_array = np.asarray(template_mass_pre,
                                         dtype='float64')
    dm.image_compose(template_array, dm.phiinvx,
                     dm.phiinvy, template_mass_pre_array)

    template_mass_pre_array *= dm.J
    W = template_mass_pre_array
    template_mass_pre = discr_space.element(W)

    # STEP 2: compute the L2 gradient
    tmpx_op = xray_trafo_op(template_mass_pre)
    tmpx_op -= noise_proj_data

    E[k+kE] += np.asarray(tmpx_op**2).sum()

    tmpx_op = xray_trafo_op.adjoint(tmpx_op)
    tmpx_array = np.array(tmpx_op, dtype='float64')
    dm.image_gradient(tmpx_array, dm.dtmpdx, dm.dtmpdy)
    dm.dtmpdx *= W
    dm.dtmpdy *= W

    dm.image_gradient(dm.sqrtJ, dm.dsqrtJdx, dm.dsqrtJdy)

    # Compute the L2 gradient of the energy functional
    density_match_L2_gradient_rec(sigma, dm.dsqrtJdx, dm.dsqrtJdy,
                                  dm.dtmpdx, dm.dtmpdy,
                                  vx, vy)

    # STEP 3:
    fftx = np.fft.fftn(vx)
    ffty = np.fft.fftn(vy)
    fftx *= dm.Linv
    ffty *= dm.Linv
    vx[:] = -np.fft.ifftn(fftx).real
    vy[:] = -np.fft.ifftn(ffty).real

    # STEP 4 (v = -grad E, so to compute the inverse
    # we solve \psiinv' = -epsilon*v o \psiinv)
    np.copyto(tmpx, vx)
    tmpx *= epsilon
    np.copyto(dm.psiinvx, dm.idx)
    dm.psiinvx -= tmpx
    # Compute forward phi also (only for output purposes)
    if dm.compute_phi:
        np.copyto(dm.psix, dm.idx)
        dm.psix += tmpx

    np.copyto(tmpy, vy)
    tmpy *= epsilon
    np.copyto(dm.psiinvy, dm.idy)
    dm.psiinvy -= tmpy
    # Compute forward phi also (only for output purposes)
    if dm.compute_phi:
        np.copyto(dm.psiy, dm.idy)
        dm.psiy += tmpy

    # STEP 5
    dm.diffeo_compose(dm.phiinvx, dm.phiinvy,
                      dm.psiinvx, dm.psiinvy,
                      tmpx, tmpy)
    np.copyto(dm.phiinvx, tmpx)
    np.copyto(dm.phiinvy, tmpy)
    # Compute forward phi also (only for output purposes)
    if dm.compute_phi:
        dm.diffeo_compose(dm.phix, dm.phiy, dm.psix, dm.psiy, tmpx, tmpy)
        np.copyto(dm.phix, tmpx)
        np.copyto(dm.phiy, tmpy)

    # STEP 6
    dm.image_compose(dm.J, dm.psiinvx, dm.psiinvy, dm.sqrtJ)
    np.copyto(dm.J, dm.sqrtJ)
    dm.divergence(vx, vy, dm.divv)
    dm.divv *= -epsilon
    np.exp(dm.divv, out=dm.sqrtJ)
    dm.J *= dm.sqrtJ
    np.sqrt(dm.J, out=dm.sqrtJ)

# Test time, set end time
end = time.clock()

# Output the computational time
print(end - start)

W = W.T
dm.J = dm.J.T
dm.phiinvx, dm.phiinvy = dm.phiinvy, dm.phiinvx
backproj = np.asarray(backproj)
backproj = backproj.T

dm.template_mass_pre = discr_space.element(dm.W.T)
rec_proj_data = xray_trafo_op(dm.template_mass_pre)

plt.figure(1, figsize=(28, 28))
plt.clf()

plt.subplot(3, 3, 1)
plt.imshow(I0, cmap='bone', vmin=dm.I0.min(), vmax=I0.max())
plt.colorbar()
plt.title('Ground truth')

plt.subplot(3, 3, 2)
plt.imshow(I1, cmap='bone', vmin=dm.I1.min(), vmax=I1.max())
plt.colorbar()
plt.title('Template')

plt.subplot(3, 3, 3)
plt.imshow(backproj, cmap='bone', vmin=backproj.min(), vmax=backproj.max())
plt.colorbar()
plt.title('Backprojection')

plt.subplot(3, 3, 4)
# plt.imshow(dm.W**2, cmap='bone', vmin=dm.I0.min(), vmax=dm.I0.max())
plt.imshow(W, cmap='bone', vmin=I1.min(), vmax=I1.max())
plt.colorbar()
plt.title('Reconstructed image by {!r} directions'.format(num_angles))
# plt.title('Warped image')

jac_ax = plt.subplot(3, 3, 5)
mycmap = 'PiYG'
# mycmap = 'Spectral'
# mycmap = 'PRGn'
# mycmap = 'BrBG'
plt.imshow(dm.J, cmap=mycmap, vmin=dm.J.min(), vmax=1.+(1.-dm.J.min()))
plt.gca().set_autoscalex_on(False)
plt.gca().set_autoscaley_on(False)
# plot_warp(dm.phiinvx, dm.phiinvy, downsample=8)
jac_colorbar = plt.colorbar()
plt.title('Jacobian')

plt.subplot(3, 3, 6)
ddmatch.plot_warp(dm.phiinvx, dm.phiinvy, downsample=4)
plt.axis('equal')
warplim = [dm.phiinvx.min(), dm.phiinvx.max(),
           dm.phiinvy.min(), dm.phiinvy.max()]
warplim[0] = min(warplim[0], warplim[2])
warplim[2] = warplim[0]
warplim[1] = max(warplim[1], warplim[3])
warplim[3] = warplim[1]

plt.axis(warplim)
# plt.axis('off')
plt.gca().invert_yaxis()
plt.gca().set_aspect('equal')
plt.title('Warp')

plt.subplot(3, 3, 7)
plt.title('stepsize = {!r}, $\sigma$ = {!r}'.format(epsilon, sigma))
plt.plot(E)
plt.ylabel('Energy')
# plt.gca().axes.yaxis.set_ticklabels(['0']+['']*8)
plt.gca().axes.yaxis.set_ticklabels([])
plt.grid(True)

plt.subplot(3, 3, 8)
plt.plot(np.asarray(proj_data)[0], 'b', np.asarray(noise_proj_data)[0], 'r')
plt.title('Theta=0, blue: truth_data, red: noisy_data, SNR = 9.17dB')
plt.gca().axes.yaxis.set_ticklabels([])
plt.axis([0, 191, -17, 32])

plt.subplot(3, 3, 9)
plt.plot(np.asarray(proj_data)[0], 'b', np.asarray(rec_proj_data)[0], 'r')
plt.title('Theta=0, blue: truth_data, red: rec result')
plt.gca().axes.yaxis.set_ticklabels([])
plt.axis([0, 191, -17, 32])
