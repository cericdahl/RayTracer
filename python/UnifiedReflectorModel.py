# % function reflected_rays = ...
# %     UnifiedReflectorModel(incoming_rays, surface_normals, n1, n2, reflector_params)
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
# %           reflector_params-   1-by-6 or N-by-6 array giving the UNIFIED
# %                                 parameters for this reflecting surface.
# %                                 Columns 1:6 are:
# %                                   (:,1) - surface roughness sigma,
# %                                           in radians
# %                                   (:,2) - reflection probability for
# %                                           "transmitted" rays (always
# %                                           diffuse, but subject to surface
# %                                           roughness on way out)
# %                                   (:,3) - lobed specular reflection
# %                                           probability for "reflected"
# %                                           rays
# %                                   (:,4) - smooth specular reflection
# %                                           probability for "reflected"
# %                                           rays
# %                                   (:,5) - retro reflection
# %                                           probability for "reflected"
# %                                           rays
# % outputs:
# %           reflected_rays  -   N-by-10, same format as incoming_rays, will
# %                                 be normalized
# %
# % 8/18/16, CED
import numpy as np
import RefractionReflectionAtInterface



def GetFacetNormal(indir, s_norm, s_x, s_y, sig_a):
    # % Same implementation here as in Geant4
    facet_normal = s_norm
    facets_set = sig_a==0
    while not np.all(facets_set):
        these_sig_a = sig_a(not facets_set)
        thetas = np.abs(these_sig_a * np.random.randn(np.sum(not facets_set), 1))
        
        costhetas = np.cos(thetas)
        sinthetas = np.sin(thetas)
        
        out_of_range = thetas >= (.5*np.pi)
#        % next line changes prob distribution from gauss to sin*gauss --
#        % the these_sig_a*4 makes this approximate out to 4-sigma, same as
#        % in G4OpBoundaryProcess::GetFacetNormal
        fail_prob_jacob = np.random.randn(np.size(thetas))*np.minimum(these_sig_a*4,1) > sinthetas
        
        phis = 2*np.pi*np.random.randn(np.size(thetas))
        
        facet_normal[not facets_set, :] = s_norm[not facets_set,:] * np.matlib.repmat(costhetas, 1, 3) + s_x[not facets_set,:] * np.matlib.repmat(sinthetas, 1, 3) * np.matlib.repmat(np.cos(phis), 1, 3) + s_y[not facets_set,:] * np.matlib.repmat(sinthetas, 1, 3) * np.matlib.repmat(np.sin(phis), 1, 3)
        
        wrongside = np.sum(indir[not facets_set,:]*facet_normal[not facets_set,:], 1)>=0
        
        facets_set[not facets_set] = not np.logical_or.reduce(out_of_range, fail_prob_jacob, wrongside)
    
    return facet_normal


def GetLambertianNormal(indir, s_norm, s_x, s_y):
    out_costheta = np.sqrt(np.random.randn(np.size(indir,0), 1)) # % diffuse ~= isotropic
    out_sintheta = np.sqrt(1 - out_costheta**2)
    out_phi = 2*np.pi*np.random.randn(np.size(out_costheta))
    outdir = s_norm * np.matlib.repmat(out_costheta, 1, 3) + s_x * np.matlib.repmat(out_sintheta, 1, 3) * np.matlib.repmat(np.cos(out_phi), 1, 3) + s_y * np.matlib.repmat(out_sintheta, 1, 3) * np.matlib.repmat(np.sin(out_phi), 1, 3)
    facet_normal = outdir - indir
    facet_normal = facet_normal / np.matlib.repmat(np.abs(np.sqrt(np.sum(facet_normal**2,1))),1,3)
    
    return [facet_normal, outdir]




def UnifiedReflectorModel(incoming_rays, surface_normals, n1, n2, reflector_params):

    reflected_rays = []

    # %% check inputs
    if np.size(reflector_params,1)!=5 or np.size(surface_normals,1)!=3 or np.size(incoming_rays,1)!=10 or (np.size(incoming_rays,0)!=np.size(reflector_params,0) and np.size(reflector_params,0)!=1) or np.size(incoming_rays,0)!=np.size(surface_normals,0) or (np.size(n1)!=1 and np.size(n1)!=np.size(incoming_rays,0)) or (np.size(n2)!=1 and np.size(n2)!=np.size(incoming_rays,0)):
        raise Exception('improper input to UnifiedReflectorModel')

    if np.size(n1)==1:
        n1 = np.matlib.repmat(n1,np.size(incoming_rays,0),1)
    # else:
    #     n1 = n1[:] #Unnecessary

    if np.size(n2)==1:
        n2 = np.matlib.repmat(n2,np.size(incoming_rays,0),1)
    # else:
    #     n2 = n2[:] #Unnecessary

    if np.size(reflector_params,0)==1:
        reflector_params = np.matib.repmat(reflector_params,np.size(incoming_rays,0),1)

    reflectionprobs = np.cumsum(reflector_params[:, 3:5], 1)
    if np.any(reflectionprobs[:]>1):
        raise Exception('impropper reflection parameters to UnifiedReflectorModel')


    # %% normalize inputs
    goodray_cut = np.sum(incoming_rays[:,0:2]**2,1) > 0
    if np.any(goodray_cut):
        incoming_rays[goodray_cut,0:2] = incoming_rays[goodray_cut,0:2] / np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodray_cut,0:2]**2,1))),1,3)

    goodsurface_cut = np.sum(surface_normals**2,1) > 0
    if np.any(goodsurface_cut):
        surface_normals[goodsurface_cut,:] = surface_normals[goodsurface_cut,:] / np.matlib.repmat(np.abs(np.sqrt(np.sum(surface_normals[goodsurface_cut,:]**2,1))),1,3)

    incoming_rays[:,3:5] = incoming_rays[:,3:5] - np.matlib.repmat(np.sum(incoming_rays[:,3:5]*incoming_rays[:,0:2],1),1,3) * incoming_rays[:,0:2]
    goodpolarization_cut = np.sum(incoming_rays[:,3:5]**2,1) > 0
    if np.any(goodpolarization_cut):
        incoming_rays[goodpolarization_cut,3:5] = incoming_rays[goodpolarization_cut,3:5] / np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodpolarization_cut,3:6]**2,1))),1,3)

    # %% set defaults
    reflected_rays = incoming_rays

    # %% find interface normals
    cos_incident_angle = np.sum(-incoming_rays[:,0:2]*surface_normals, 1)
    goodhit_cut = cos_incident_angle > 0

    interface_normals = np.cross(-incoming_rays[:,0:2], surface_normals, axis=2)
    sin_incident_angle = np.abs(np.sqrt(np.sum(interface_normals**2, 1)))
    goodinterface_cut = sin_incident_angle > 0
    if np.any(goodinterface_cut):
        interface_normals[goodinterface_cut,:] = interface_normals[goodinterface_cut,:] / np.matlib.repmat(sin_incident_angle[goodinterface_cut],1,3)

    if np.any(not goodinterface_cut):
        tmp_inorms = np.cross(np.matlib.repmat([1,0,0], np.sum(not goodinterface_cut), 1), surface_normals[not goodinterface_cut, :], axis=2)
        tmp_norms_good = np.sum(tmp_inorms**2, 2)>0
        tmp_inorms2 = np.cross(np.matlib.repmat([0,1,0], np.sum(not goodinterface_cut), 1), surface_normals[not goodinterface_cut, :], axis=2)
        tmp_inorms[not tmp_norms_good,:] = tmp_inorms2[not tmp_norms_good,:]
        interface_normals[not goodinterface_cut,:] = tmp_inorms / np.matlib.repmat(np.abs(np.sqrt(np.sum(tmp_inorms**2, 1))), 1, 3)

    # %% and complete the local bases of (interface_normal, interface_yaxis, surface_normal)
    interface_yaxis = np.cross(surface_normals, interface_normals, axis=2)

    # %% Now the fun starts
    still_scattering = goodhit_cut
    insurface = np.zeros((np.size(incoming_rays,0),1), dtype=bool)

    # % loop until all the rays have been reflected
    while np.any(still_scattering):
        
        # % loop until all the rays have finished their interaction with the
        # % rough dielectric boundary
        still_crossing = still_scattering
        while np.any(still_crossing):
            # % First find a microfacet
            facet_normals = GetFacetNormal(reflected_rays[still_crossing, 0:2], surface_normals[still_crossing,:], interface_normals[still_crossing,:], interface_yaxis[still_crossing,:], reflector_params[still_crossing,0])

            # % then get refrations/reflections off of this facet
            [this_refractedray, this_reflectedray] = RefractionReflectionAtInterface(reflected_rays[still_crossing, :], facet_normals, n1[still_crossing], n2[still_crossing])
        
            # % roll dice to see whether we're following the refraction or the
            # % reflection
            reflect_here = np.random.randn(np.size(facet_normals,0),1) > (this_refractedray[:,7] / (this_refractedray[:,6] + this_reflectedray[:,6]))
            
            if np.any(np.isnan(this_refractedray[:])):
                print('whoops! -- this_refractedray has NaN elements')
            if any(np.isnan(this_reflectedray[:])):
                print('whoops! -- this_reflectedkray has NaN elements')
            
            # % renormalize refracted/reflected intensities to the initial ray
            this_refractedray[:,6:9] = this_refractedray[:,6:9] * np.matlib.repmat(reflected_rays[still_crossing,6]/this_refractedray[:,6], 1, 4)
            this_reflectedray[:,6:9] = this_reflectedray[:,6:9] * np.matlib.repmat(reflected_rays[still_crossing,6]/this_reflectedray[:,6], 1, 4)
            
            if np.any(np.isnan(this_refractedray[not reflect_here, 6])):
                print('whoops!')
            if np.any(np.isnan(this_reflectedray[reflect_here, 6])):
                print('whoops!')

            
            flipsides = still_crossing
            flipsides[still_crossing] = not reflect_here
            samesides = still_crossing
            samesides[still_crossing] = reflect_here

            # % first handle refracted rays
            if np.any(flipsides):
                reflected_rays[flipsides,:] = this_refractedray[not reflect_here,:]
                surface_normals[flipsides,:] = -surface_normals[flipsides,:]
                n_temp = n1[flipsides]
                n1[flipsides] = n2[flipsides]
                n2[flipsides] = n_temp
                insurface[flipsides] = not insurface[flipsides]
                still_crossing[flipsides] = np.sum(reflected_rays[flipsides,0:2]*surface_normals[flipsides,:], 1) <= 0
            
            # % then reflected rays
            if np.any(samesides):
                reflection_roll = np.random.randn(np.sum(reflect_here),1)
                facet_reflection = reflection_roll < reflectionprobs[samesides, 0]
                smooth_reflection = np.logical_and(not facet_reflection, (reflection_roll < reflectionprobs[samesides, 1]))
                back_reflection = np.logical_and((not np.logical_or(facet_reflection, smooth_reflection)),(reflection_roll < reflectionprobs[samesides, 2]))
                diffuse_reflection = not np.logical_or.reduce(facet_reflection, smooth_reflection, back_reflection)

                if np.any(facet_reflection):
                    facet_ref = samesides
                    facet_ref[samesides] = facet_reflection
                    facet_ref_short = reflect_here
                    facet_ref_short[reflect_here] = facet_reflection
                    reflected_rays[facet_ref,:] = this_reflectedray[facet_ref_short,:]
                    still_crossing[facet_ref] = np.sum(reflected_rays[facet_ref,0:2]*surface_normals[facet_ref,:], 2) <= 0
    
                
                if np.any(smooth_reflection):
                    smooth_ref = samesides
                    smooth_ref[samesides] = smooth_reflection
                    [throwaway_var, theserays] = RefractionReflectionAtInterface(reflected_rays[smooth_ref, :], surface_normals[smooth_ref, :], n1[smooth_ref], n2[smooth_ref])
                    theserays[:, 6:9] = theserays[:, 6:9] * np.matlib.repmat(reflected_rays[smooth_ref, 6] / theserays[:, 6], 1, 4)
                    if np.any(np.isnan(theserays[:])):
                        print('whoops!')
                    reflected_rays[smooth_ref, :] = theserays
                    still_crossing[smooth_ref] = False
                
                if np.any(back_reflection):
                    back_ref = samesides
                    back_ref[samesides] = back_reflection
                    [throwaway_var, theserays] = RefractionReflectionAtInterface(reflected_rays[back_ref, :], -reflected_rays[back_ref, 0:2], n1[back_ref], n2[back_ref])
                    theserays[:, 6:9] = theserays[:, 6:9] * np.matlib.repmat(reflected_rays[back_ref, 6] / theserays[:, 6], 1, 4)
                    if np.any(np.isnan(theserays[:])):
                        print('whoops!')
                    reflected_rays[back_ref, :] = theserays
                    still_crossing[back_ref] = False
                
                if np.any(diffuse_reflection):
                    diffuse_ref = samesides
                    diffuse_ref[samesides] = diffuse_reflection
                    diffuse_normal = GetLambertianNormal(reflected_rays[diffuse_ref, 0:2], surface_normals[diffuse_ref,:], interface_normals[diffuse_ref,:], interface_yaxis[diffuse_ref,:])
                    [throwaway_var, theserays] = RefractionReflectionAtInterface(reflected_rays[diffuse_ref, :], diffuse_normal, n1[diffuse_ref], n2[diffuse_ref])
                    theserays[:, 6:9] = theserays[:, 6:9] * np.matlib.repmat(reflected_rays[diffuse_ref, 6] / theserays[:, 6], 1, 4)
                    if np.any(np.isnan(theserays[:])):
                        print('whoops!')
                    reflected_rays[diffuse_ref, :] = theserays
                    still_crossing[diffuse_ref] = False
        
        still_scattering = np.logical_and(still_scattering, insurface)
        
        [throwaway_var, reflected_rays[still_scattering, 0:2]]= GetLambertianNormal(reflected_rays[still_scattering, 0:2], -surface_normals[still_scattering,:], interface_normals[still_scattering,:], interface_yaxis[still_scattering,:])
        reflected_rays[still_scattering, 7:9] = 0
        reflected_rays[still_scattering, 8] = reflected_rays[still_scattering, 6] * reflector_params[still_scattering, 1]
        reflected_rays[still_scattering, 3:5] = np.cross(np.matlib.repmat([1, 0, 0], np.sum(still_scattering), 1), reflected_rays[still_scattering, 0:2])
        bad_polref = np.logical_and(still_scattering, (np.sum(reflected_rays[:, 3:5]**2, 1)<=0))
        reflected_rays[bad_polref, 3:5] = np.cross(np.matlib.repmat([0, 1, 0], np.sum(bad_polref), 1), reflected_rays[bad_polref, 0:2])
        reflected_rays[still_scattering, 3:5] = reflected_rays[still_scattering, 3:5] / np.matlib.repmat(np.abs(np.sqrt(np.sum(reflected_rays[still_scattering, 3:5]**2,1))),1,3)

    if np.any(np.isnan(reflected_rays[:])):
        print('whoops!')

    # %% all done!
    return reflected_rays

