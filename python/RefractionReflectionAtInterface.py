# % function [refracted_rays reflected_rays] = ...
# %     RefractionReflectionAtInterface(incoming_rays, surface_normals, n1, n2, tir_handling)
# %
# % inputs:
# %           incoming_rays   -   N-by-10 vector, where N is the number of
# %                                 rays to refract.  The first three
# %                                 elements are the x, y, and z components
# %                                 of the forward direction of the ray (will
# %                                 be normalized if not already), the 4-6
# %                                 elements define the polarization
# %                                 reference frame (will be projected into
# %                                 plane perpendicular to ray and normalized
# %                                 if not already), and the 7-10 elements
# %                                 are the stokes parameters for the ray (7
# %                                 gives intensity, 8-10 give polarization)
# %           surface_normals -   N-by-3 vector giving the outward-pointing
# %                                 normals on the incident surface
# %           n1              -   1 or N element vector giving the index of
# %                                 refraction of the medium the incoming ray
# %                                 travels through
# %           n2              -   1 or N element vector giving the index of
# %                                 refraction of the medium the on the other
# %                                 side of the interface
# %           tir_handling    -   This determines what the refracted_rays
# %                                 output is in the case of total internal
# %                                 reflection.  The default (-1) gives a
# %                                 refracted ray tangent to the surface with
# %                                 zero intensity.  Any value >=0 will give
# %                                 a ray with the same direction and
# %                                 polarization as the reflected ray, with
# %                                 intensity equal to the reflected
# %                                 intensity times tir_handling.  This lets
# %                                 you treat tir-rays like refracted rays,
# %                                 which can be handy in geometry sims.
# % outputs:
# %           refracted_rays  -   N-by-10, same format as incoming_rays, will
# %                                 be normalized
# %           reflected_rays  -   N-by-10, same format as incoming_rays, will
# %                                 be normalized
# %
# % 12/15/09, CED

import numpy as np
import math

def RefractionReflectionAtInterface(incoming_rays, surface_normals, n1, n2, tir_handling = -1):

    refracted_rays = []
    reflected_rays = []

    # %% check inputs
    if np.size(tir_handling) == 0:
        tir_handling = -1

    if np.size(surface_normals,1)!=3 or np.size(incoming_rays,1)!=10 or np.size(incoming_rays,0)!=np.size(surface_normals,0) or (np.size(n1)!=1 and np.size(n1)!=np.size(incoming_rays,0)) or (np.size(n2)!=1 and np.size(n2)!=np.size(incoming_rays,0)) or (np.size(tir_handling)!=1 and np.size(tir_handling)!=np.size(incoming_rays,0)):
        raise Exception('impropper input to RefractionReflactionAtInterface')
        

    if np.size(n1)==1:
        n1 = np.matlib.repmat(n1,np.size(incoming_rays,0),1)
    # else:
    #     n1 = n1[:]  #Unnecessary

    if np.size(n2)==1:
        n2 = np.matlib.repmat(n2,np.size(incoming_rays,0),1)
    # else:
    #     n2 = n2[:]  #Unnecessary


    if np.size(tir_handling)==1:
        tir_handling = np.matlib.repmat(tir_handling,np.size(incoming_rays,0),1)
    # else:
    #     tir_handling = tir_handling[:]   #Unnecessary


    # %% normalize inputs
    goodray_cut = np.sum(incoming_rays[:,0:3]**2,1) > 0
    if np.any(goodray_cut):
        incoming_rays[goodray_cut,0:3] = incoming_rays[goodray_cut,0:3] / np.swapaxes(np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodray_cut,0:3]**2,1))),3,1), 1, 0) # swapped 3 and 1, changed index [0:2] --> [0:3], swapped axes (3,1000) --> (1000,3)


    goodsurface_cut = np.sum(surface_normals**2,1) > 0
    if np.any(goodsurface_cut):
        surface_normals[goodsurface_cut,:] = surface_normals[goodsurface_cut,:] / np.swapaxes(np.matlib.repmat(np.abs(np.sqrt(np.sum(surface_normals[goodsurface_cut,:]**2,1))),3,1), 1, 0)

    incoming_rays[:,3:6] = incoming_rays[:,3:6] - np.swapaxes(np.matlib.repmat(np.sum(incoming_rays[:,3:6]*incoming_rays[:,0:3],1),3,1), 1, 0) * incoming_rays[:,0:3]
    goodpolarization_cut = np.sum(incoming_rays[:,3:6]**2,1) > 0
    if np.any(goodpolarization_cut):
        incoming_rays[goodpolarization_cut,3:6] = incoming_rays[goodpolarization_cut,3:6] / np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodpolarization_cut,3:6]**2,1))),1,3)
    #print(incoming_rays[goodpolarization_cut,3:6].shape)
    #print(np.abs(np.sqrt(np.sum(incoming_rays[goodpolarization_cut,3:6]**2,1))).shape)
    # %% set defaults
    refracted_rays = incoming_rays
    reflected_rays = incoming_rays
    reflected_rays[:,0:2] = -reflected_rays[:,0:2]
    reflected_rays[:,6:9] = 0

    # %% find interface normals
    cos_incident_angle = np.sum(-incoming_rays[:,0:3]*surface_normals, 1)
    goodhit_cut = cos_incident_angle > 0

    interface_normals = np.cross(-incoming_rays[:,0:3], surface_normals, axis=1)
    sin_incident_angle = np.abs(np.sqrt(np.sum(interface_normals**2, 1)))
    goodinterface_cut = sin_incident_angle > 0
    if np.any(goodinterface_cut):
        interface_normals[goodinterface_cut,:] = interface_normals[goodinterface_cut,:] / np.swapaxes(np.matlib.repmat(sin_incident_angle[goodinterface_cut],3,1), 1, 0)

    # %% rotate stokes parameters to new basis
    c_rot = np.sum(interface_normals * incoming_rays[:,3:6], 1)
    s_rot = np.sum(np.cross(interface_normals,incoming_rays[:,3:6],axis=1) * incoming_rays[:,0:3], 1)
    c2_rot = c_rot**2 - s_rot**2
    s2_rot = 2*c_rot*s_rot
    # % for scalar c2_rot and s2_rot, we would write:
    # %     rotmat = [c2_rot s2_rot ; -s2_rot c2_rot];
    # % but the transformation is different for every ray, so.......
    if np.any(goodinterface_cut):
        incoming_rays[goodinterface_cut,3:6] = interface_normals[goodinterface_cut,:]
    # % again, if we had a 2x2 rotmat, we would write:
    # %     incoming_rays(goodinterface_cut,8:9) = incoming_rays(goodinterface_cut,8:9) * rotmat;
    # % but instead we do it one element at a time
        old_polarization = incoming_rays[goodinterface_cut,7:9]
        incoming_rays[goodinterface_cut,7] = old_polarization[:,0]*c2_rot[goodinterface_cut] - old_polarization[:,1]*s2_rot[goodinterface_cut]
        incoming_rays[goodinterface_cut,8] = old_polarization[:,0]*s2_rot[goodinterface_cut] + old_polarization[:,1]*c2_rot[goodinterface_cut]
        refracted_rays[goodinterface_cut,3:9] = incoming_rays[goodinterface_cut,3:9]
        reflected_rays[goodinterface_cut,3:5] = incoming_rays[goodinterface_cut,3:5]

    # %% find complex amplitudes for incoming rays in new basis
    # % We break the incoming ray into 3 polarized planes waves -- one in the
    # % polarization direction indicated by the stokes parameters, one in the
    # % interface_normal, and one in the interface-plane -- the latter two make
    # % up the unpolarized portion of the incoming ray
    amplitudes = np.zeros((incoming_rays.shape[0],3,2))
    print(1j*np.arctan2(incoming_rays[:,9],incoming_rays[:,8]))
    polarized_intensity = np.abs(np.sqrt(np.sum(incoming_rays[:,7:9]**2,1)))
    amplitudes[:,0,0] = np.abs(np.sqrt(.5 * (polarized_intensity + incoming_rays[:,7])))
    amplitudes[:,0,1] = np.abs(np.sqrt(.5 * (polarized_intensity - incoming_rays[:,7]))) * np.exp(1j*math.atan2(incoming_rays[:,9],incoming_rays[:,8]))
    amplitudes[:,1,0] = np.abs(np.sqrt(.5 * (incoming_rays[:,6] - polarized_intensity)))
    amplitudes[:,2,1] = amplitudes[:,1,0]

    # %% calculated reflected and refracted amplitudes
    sin_refracted_angle = sin_incident_angle * n1 / n2
    cos_refracted_angle = np.sqrt(1 - sin_refracted_angle**2)

    rs = (n1*cos_incident_angle - n2*cos_refracted_angle) / (n1*cos_incident_angle + n2*cos_refracted_angle)
    rp = -(n1*cos_refracted_angle - n2*cos_incident_angle) / (n1*cos_refracted_angle + n2*cos_incident_angle)

    rs[np.logical_or(n2==np.nf, n2==-np.inf)] = -1
    rp[np.logical_or(n2==np.inf, n2==-np.inf)] = 1

    ts = np.abs(np.sqrt(1-np.conj(rs)*rs))
    tp = np.abs(np.sqrt(1-np.conj(rp)*rp))
    refracted_amplitudes = amplitudes * np.tile(np.reshape([ts, tp],(-1,1,2)),(1, 3, 1))
    reflected_amplitudes = amplitudes * np.tile(np.reshape([rs, rp],(-1,1,2)),(1, 3, 1))

    # %% get back to stokes parameters
    if np.any(goodhit_cut):
        refracted_rays[goodhit_cut,6] = np.sum(np.sum(np.conj(refracted_amplitudes[goodhit_cut,:,:]) * refracted_amplitudes[goodhit_cut,:,:],2),1)
        refracted_rays[goodhit_cut,8] = -np.sum(np.diff(np.conj(refracted_amplitudes[goodhit_cut,:,:]) * refracted_amplitudes[goodhit_cut,:,:],1,3),1)
        refracted_rays[goodhit_cut,9] = np.sum(np.real(2 * np.conj(refracted_amplitudes[goodhit_cut,:,0]) * refracted_amplitudes[goodhit_cut,:,1]),1)
        refracted_rays[goodhit_cut,10] = np.sum(np.imag(2 * np.conj(refracted_amplitudes[goodhit_cut,:,0]) * refracted_amplitudes[goodhit_cut,:,1]),1)

        reflected_rays[goodhit_cut,6] = np.sum(np.sum(np.conj(reflected_amplitudes[goodhit_cut,:,:]) * reflected_amplitudes[goodhit_cut,:,:],2),1)
        reflected_rays[goodhit_cut,7] = -np.sum(np.diff(np.conj(reflected_amplitudes[goodhit_cut,:,:]) * reflected_amplitudes[goodhit_cut,:,:],1,3),1)
        reflected_rays[goodhit_cut,8] = np.sum(np.real(2 * np.conj(reflected_amplitudes[goodhit_cut,:,0]) * reflected_amplitudes[goodhit_cut,:,1]),1)
        reflected_rays[goodhit_cut,9] = np.sum(np.imag(2 * np.conj(reflected_amplitudes[goodhit_cut,:,0]) * reflected_amplitudes[goodhit_cut,:,1]),1)

    # %% consider surface_normal to be -x, interface_normal to be +z
    new_yaxis = np.cross(surface_normals, interface_normals, axis=2)
    goodcut = np.logical_and(goodhit_cut, np.sum(new_yaxis**2,1) > 0)

    if np.any(goodcut):
        new_yaxis[goodcut,:] = new_yaxis[goodcut,:] / np.matlib.repmat(np.abs(np.sqrt(np.sum(new_yaxis[goodcut,:]**2,1))),1,3)

    # %% calculate reflected ray direction for non-normal incidence
    if np.any(goodcut):
        reflected_rays[goodcut,0:2] = np.matlib.repmat(cos_incident_angle[goodcut],1,3) * surface_normals[goodcut,:] - np.matlib.repmat(sin_incident_angle[goodcut],1,3) * new_yaxis[goodcut,:]

    # %% calculated refracted ray direction
    total_internal_reflection_cut = np.logical_and(goodcut, sin_refracted_angle) >= 1
    refracted_cut = np.logical_and(goodcut, not total_internal_reflection_cut)

    if np.any(refracted_cut):
        refracted_rays[refracted_cut, 0:2] = -np.matlib.repmat(cos_refracted_angle[refracted_cut],1,3) * surface_normals[refracted_cut,:] - np.matlib.repmat(sin_refracted_angle[refracted_cut],1,3) * new_yaxis[refracted_cut,:]

    if np.any(total_internal_reflection_cut):
        if np.any(np.logical_and(total_internal_reflection_cut, (tir_handling<0))):
            refracted_rays[np.logical_and(total_internal_reflection_cut, (tir_handling<0)), 0:2] = -new_yaxis[np.logical_and(total_internal_reflection_cut, (tir_handling<0)),:]

        if np.any(np.logical_and(total_internal_reflection_cut, (tir_handling>=0))):
            refracted_rays[np.logical_and(total_internal_reflection_cut, (tir_handling>=0)), 0:5] = reflected_rays[np.logical_and(total_internal_reflection_cut, (tir_handling>=0)),0:5]
            refracted_rays[np.logical_and(total_internal_reflection_cut, (tir_handling>=0)), 6:9] = reflected_rays[np.logical_and(total_internal_reflection_cut, (tir_handling>=0)),6:9] * np.matlib.repmat(tir_handling[np.logical_and(total_internal_reflection_cut, (tir_handling>=0))],1,4)

    # %% all done!
    return [refracted_rays, reflected_rays]
    
