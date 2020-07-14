#% function scattered_rays = ...
#%     RayleighScattering(incoming_rays)
#%
#% This function Rayleigh scatters each of the incoming rays.  The output
#% direction is randomized according to the Rayleigh scattering distribution
#% (dependent on the polarization of each incoming ray) and the output
#% polarization (i.e. Stokes parameters) is fixed based on the incoming
#% polarization and scattering direction.  Intensity of the ray is
#% unchanged.  For a nice description of Rayleigh scattering in terms of
#% Stokes parameters, see https://www.astro.umd.edu/~jph/notes3.pdf
#%
#% inputs:
#%           incoming_rays   -   N-by-10 vector, where N is the number of
#%                                 rays to refract.  The first three
#%                                 elements are the x, y, and z components
#%                                 of the forward direction of the ray (will
#%                                 be normalized if not already), the 4-6
#%                                 elements define the polarization
#%                                 reference frame (will be projected into
#%                                 plane perpendicular to ray and normalized
#%                                 if not already), and the 7-10 elements
#%                                 are the stokes parameters for the ray (7
#%                                 gives intensity, 8-10 give polarization)
#%           azimuth_precision - The plane of scattering is determined by
#%                                 a dice roll to a lookup table of length
#%                                 (azimuth_precision) with equally-spaced
#%                                 entries in cumulative probability (and
#%                                 linear interpolation between entries).
#%                                 Default value is 1e5.
#% outputs:
#%           scattered_rays  -   N-by-10, same format as incoming_rays, will
#%                                 be normalized
#%
#%
#% 8/18/16, CED

import numpy as np
import math

#function scattered_rays = RayleighScattering(incoming_rays, azimuth_precision)
class RayleighScatteringClass(object):

    #define persistent rayleigh_azimuth
    def __init__(self):
        self.rayleigh_azimuth = []
        
    def RayleighScattering(self, incoming_rays=[], azimuth_precision = []):
#        %% set defaults
        scattered_rays = []

        if np.size(azimuth_precision) == 0:
            azimuth_precision = 1e5

#        %% check inputs
        if np.size(incoming_rays,1)!=10:
            raise Exception('impropper input to RayleighScattering')
            
#        %% initialize self.rayleigh_azimuth
        if np.size(self.rayleigh_azimuth)==0 or np.size(self.rayleigh_azimuth) != azimuth_precision:
            cumd = np.linspace(0, 1, azimuth_precision+1)
            phid = np.linspace(0, 2*np.pi, azimuth_precision*10)
            cumd_phi = (phid - .25*np.sin(2*phid)) / (2*np.pi)
            self.rayleigh_azimuth = np.transpose(np.interp(cumd_phi, cumd,phid)) #switched 2nd and 2rd inputs for matlab interp1 -> python interp
            
#        %% normalize inputs
        goodray_cut = np.sum(incoming_rays[:,0:2]**2,1) > 0
        if np.any(goodray_cut):
            incoming_rays[goodray_cut,0:2] = incoming_rays[goodray_cut,0:2] / np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodray_cut,0:2]**2,1))),1,3)

        incoming_rays[:,3:5] = incoming_rays[:,3:5] - np.matlib.repmat(np.sum(incoming_rays[:,3:5]*incoming_rays[:,0:2],1),1,3) * incoming_rays[:,0:2]
        goodpolarization_cut = np.sum(incoming_rays[:,3:5]**2,1) > 0
        if np.any(goodpolarization_cut):
            incoming_rays[goodpolarization_cut,3:5] = incoming_rays[goodpolarization_cut,3:5] / np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodpolarization_cut,3:5]**2,1))),1,3)
        
        

#        %% roll dice to determine scattering angle
        scatter_dice = np.random.randn(np.size(incoming_rays,0), 3)
#        % column 1 for the lin-polarized vs not-lin-polarized bit
#        % column 2 for phi
#        % column 3 for theta

#        %% Determine scattering plane relative to polarization reference
        linpol_intensity = np.abs(np.sqrt(np.sum(incoming_rays[:,7:8]**2,1)))
        linpol_frac = linpol_intensity / incoming_rays[:,6]
        linpol_scatter = linpol_frac > scatter_dice[:, 0]

        scatter_phi = 2*np.pi*scatter_dice[:, 1]

        if np.any(linpol_scatter):
            linpol_angle = .5 * math.atan2(incoming_rays[linpol_scatter,8], incoming_rays[linpol_scatter,7])
            ix = math.floor(azimuth_precision*scatter_dice[linpol_scatter, 2]) + 1
            ix = np.minimum(ix, azimuth_precision)
            ix_frac = azimuth_precision*scatter_dice[linpol_scatter, 1] - (ix-1)
            scatter_phi[linpol_scatter] = linpol_angle + self.rayleigh_azimuth[ix] + ix_frac*(self.rayleigh_azimuth[ix+1] - self.rayleigh_azimuth[ix])

#        %% Now rotate stokes parameters to the scattering plane
        c_rot = np.cos(scatter_phi)
        s_rot = np.sin(scatter_phi)
        c2_rot = c_rot**2 - s_rot**2
        s2_rot = 2 * c_rot * s_rot

        old_polarization = incoming_rays[:,7:8]
        incoming_rays[:,7] = old_polarization[:,0] * c2_rot + old_polarization[:,1] * s2_rot
        incoming_rays[:,8] = -old_polarization[:,1] * s2_rot + old_polarization[:,2] * c2_rot

        old_reference = incoming_rays[:,3:5]
        old_reference_perp = np.cross(incoming_rays[:, 1:2], incoming_rays[:, 3:5])
        incoming_rays[:, 3:5] = old_reference * np.matlib.repmat(c_rot, 1, 3) + old_reference_perp * np.matlib.repmat(s_rot, 1, 3)

        new_reference_perp = np.cross(incoming_rays[:, 0:2], incoming_rays[:, 3:5])

#        %% Now determine scattering angle
        a = incoming_rays[:, 7] / incoming_rays[:, 6] # % -1 < a < 1
#        % now cos_theta is the solution of a cubic polynomial with coefficients
#        % determined by a and scatter_dice(:, 3)
        z = (2-a) * (1 - 2*scatter_dice[:, 2])
        s = np.sqrt(z**2 + (1-a)**3/(1+a))
        A = np.sign(z+s) * np.abs(z + s)**(1/3)
        B = np.sign(z-s) * np.abs(z - s)**(1/3)
        cos_theta = (A + B) * (1 + a)**(-1/3)
        sin_theta = np.sqrt(1 - cos_theta**2)

#        %% now determine new direction
        scattered_rays = incoming_rays
        scattered_rays[:,0:2] = incoming_rays[:,0:2] * np.matlib.repmat(cos_theta, 1, 3) + old_reference * np.matlib.repmat(sin_theta, 1, 3) * np.matlib.repmat(c_rot, 1, 3) + old_reference_perp * np.matlib.repmat(sin_theta, 1, 3) * np.matlib.repmat(s_rot, 1, 3)
        scattered_rays[:, 3:5] = np.cross(new_reference_perp, scattered_rays[:, 0:2])

#        %% and new ray polarization
        c2theta = cos_theta**2

        R11 = c2theta + 1
        R12 = c2theta - 1

        rnorm = incoming_rays[:, 6] / (incoming_rays[:, 6] * R11 + incoming_rays[:, 7] * R12)
        scattered_rays[:, 7] = (incoming_rays[:, 6] * R12 + incoming_rays[:, 7] * R11) * rnorm

        R33 = 2*cos_theta * rnorm

        scattered_rays[:, 8] = incoming_rays[:, 8] * R33
        scattered_rays[:, 9] = incoming_rays[:, 9] * R33

        return scattered_rays
