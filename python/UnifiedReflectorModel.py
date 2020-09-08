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
    facets_set = sig_a == 0
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
    out_costheta = np.sqrt(np.random.random_sample((np.size(indir,0), 1))) # % diffuse ~= isotropic | randn --> random_sample for uniform distribution
    out_sintheta = np.sqrt(1 - out_costheta**2)
    out_phi = 2*np.pi*np.random.random_sample(out_costheta.shape) # randn --> random_sample for uniform distribution
    outdir = s_norm * np.matlib.repmat(out_costheta, 1, 3) + s_x * np.matlib.repmat(out_sintheta, 1, 3) * np.matlib.repmat(np.cos(out_phi), 1, 3) + s_y * np.matlib.repmat(out_sintheta, 1, 3) * np.matlib.repmat(np.sin(out_phi), 1, 3)
    facet_normal = outdir - indir
    facet_normal = facet_normal / np.matlib.repmat(np.abs(np.sqrt(np.sum(facet_normal**2,1)))[:,np.newaxis],1,3)
    
    return [facet_normal, outdir]




def UnifiedReflectorModel(incoming_rays, surface_normals, n1, n2, reflector_params):

    reflected_rays = []

    # %% check inputs
    if np.size(reflector_params,1)!=5 or np.size(surface_normals,1)!=3 or np.size(incoming_rays,1)!=10 or (np.size(incoming_rays,0)!=np.size(reflector_params,0) and np.size(reflector_params,0)!=1) or np.size(incoming_rays,0)!=np.size(surface_normals,0) or (np.size(n1)!=1 and np.size(n1)!=np.size(incoming_rays,0)) or (np.size(n2)!=1 and np.size(n2)!=np.size(incoming_rays,0)):
        raise Exception('improper input to UnifiedReflectorModel')

    if np.size(n1)==1:
        n1 = np.matlib.repmat(n1[:,np.newaxis],np.size(incoming_rays,0),1) # is the np.newaxis necessary here?
    # else:
    #     n1 = n1[:] #Unnecessary

    if np.size(n2)==1:
        n2 = np.matlib.repmat(n2,np.size(incoming_rays,0),1)
    # else:
    #     n2 = n2[:] #Unnecessary

    if np.size(reflector_params,0)==1:
        reflector_params = np.matlib.repmat(reflector_params,np.size(incoming_rays,0),1)

    reflectionprobs = np.cumsum(reflector_params[:, 2:5], axis=1)
    if np.any(reflectionprobs[:]>1):
        raise Exception('impropper reflection parameters to UnifiedReflectorModel')


    # %% normalize inputs
    goodray_cut = np.sum(incoming_rays[:,0:3]**2,1) > 0
    if np.any(goodray_cut):
        incoming_rays[goodray_cut,0:3] = incoming_rays[goodray_cut,0:3] / np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodray_cut,0:3]**2,1)))[:, np.newaxis],1,3)

    goodsurface_cut = np.sum(surface_normals**2,1) > 0
    if np.any(goodsurface_cut):
        surface_normals[goodsurface_cut,:] = surface_normals[goodsurface_cut,:] / np.matlib.repmat(np.abs(np.sqrt(np.sum(surface_normals[goodsurface_cut,:]**2,1)))[:, np.newaxis],1,3)

    incoming_rays[:,3:6] = incoming_rays[:,3:6] - np.matlib.repmat(np.sum(incoming_rays[:,3:6]*incoming_rays[:,0:3],1)[:, np.newaxis],1,3) * incoming_rays[:,0:3]
    goodpolarization_cut = np.sum(incoming_rays[:,3:6]**2, axis=1) > 0
    if np.any(goodpolarization_cut):
        incoming_rays[goodpolarization_cut,3:6] = incoming_rays[goodpolarization_cut,3:6] / np.matlib.repmat(np.abs(np.sqrt(np.sum(incoming_rays[goodpolarization_cut,3:6]**2,1)))[:, np.newaxis],1,3)

    # %% set defaults
    reflected_rays = np.copy(incoming_rays)
    #print("reflected check 1: " + str(reflected_rays))

    # %% find interface normals
    cos_incident_angle = np.sum(-incoming_rays[:,0:3] * surface_normals, axis=1)
    goodhit_cut = cos_incident_angle > 0

    interface_normals = np.cross(-incoming_rays[:,0:3], surface_normals, axis=1)
    sin_incident_angle = np.abs(np.sqrt(np.sum(interface_normals**2, 1)))
    goodinterface_cut = sin_incident_angle > 0
    if np.any(goodinterface_cut):
        interface_normals[goodinterface_cut,:] = interface_normals[goodinterface_cut,:] / np.matlib.repmat(sin_incident_angle[goodinterface_cut, np.newaxis],1,3)

    if np.any(~goodinterface_cut):
        tmp_inorms = np.cross(np.matlib.repmat([1,0,0], np.sum(not goodinterface_cut), 1), surface_normals[not goodinterface_cut, :], axis=1)
        tmp_norms_good = np.sum(tmp_inorms**2, 2)>0
        tmp_inorms2 = np.cross(np.matlib.repmat([0,1,0], np.sum(not goodinterface_cut), 1), surface_normals[not goodinterface_cut, :], axis=1)
        tmp_inorms[not tmp_norms_good,:] = tmp_inorms2[not tmp_norms_good,:]
        interface_normals[not goodinterface_cut,:] = tmp_inorms / np.matlib.repmat(np.abs(np.sqrt(np.sum(tmp_inorms**2, 1))), 1, 3)

    # %% and complete the local bases of (interface_normal, interface_yaxis, surface_normal)
    interface_yaxis = np.cross(surface_normals, interface_normals, axis=1)

    # %% Now the fun starts
    still_scattering = np.copy(goodhit_cut)
    insurface = np.zeros((np.size(incoming_rays,0),1), dtype=bool)

    # % loop until all the rays have been reflected
    while np.any(still_scattering):
        
        # % loop until all the rays have finished their interaction with the
        # % rough dielectric boundary
        still_crossing = np.copy(still_scattering)
        while np.any(still_crossing):
            # % First find a microfacet
            print("reflected check 2: " + str(reflected_rays))
            facet_normals = GetFacetNormal(reflected_rays[still_crossing, 0:3], surface_normals[still_crossing,:], interface_normals[still_crossing,:], interface_yaxis[still_crossing,:], reflector_params[still_crossing,0])

            #print(np.any(np.isnan(reflected_rays[still_crossing, :])))
            # % then get refrations/reflections off of this facet   | getting NaNs from here, reflected_rays[still_crossing, :] has Nan(s) after first iteration
            [this_refractedray, this_reflectedray] = RefractionReflectionAtInterface.RefractionReflectionAtInterface(reflected_rays[still_crossing, :], facet_normals, n1[still_crossing], n2[still_crossing])
            print("this refracted: " + str(this_refractedray))
            print("this reflected: " + str(this_reflectedray))
        
            # % roll dice to see whether we're following the refraction or the
            # % reflection
            reflect_here = (np.random.random_sample((np.size(facet_normals,0),1)) > (this_refractedray[:,6] / (this_refractedray[:,6] + this_reflectedray[:,6]))[:,np.newaxis]).flatten() # for indexing
            
            if np.any(np.isnan(this_refractedray[:])):
                print('whoops! -- this_refractedray has NaN elements')
            if np.any(np.isnan(this_reflectedray[:])):
                print('whoops! -- this_reflectedkray has NaN elements')
            
            # % renormalize refracted/reflected intensities to the initial ray  | causing NaNs; 0 * inf = nan | nan_to_num?
            this_refractedray[:,6:10] = this_refractedray[:,6:10] * np.matlib.repmat((reflected_rays[still_crossing,6]/this_refractedray[:,6])[:,np.newaxis], 1, 4)
            print("this refracted norm: " + str(this_refractedray))
            print("divide check: " + str(np.matlib.repmat(np.divide(reflected_rays[still_crossing,6], this_reflectedray[:,6])[:,np.newaxis], 1, 4)))
            this_reflectedray[:,6:10] = this_reflectedray[:,6:10] * np.matlib.repmat((reflected_rays[still_crossing,6]/this_reflectedray[:,6])[:,np.newaxis], 1, 4) # here, s0 this_reflected=0
            print("this reflected norm: " + str(this_reflectedray))
            
            if np.any(np.isnan(this_refractedray[~reflect_here, 6])): # this is being called
                print('whoops refracted!')
            if np.any(np.isnan(this_reflectedray[reflect_here, 6])):
                print('whoops reflected!')

            
            flipsides = np.copy(still_crossing)
            flipsides[still_crossing] = ~reflect_here
            samesides = np.copy(still_crossing)
            samesides[still_crossing] = reflect_here

            # % first handle refracted rays
            if np.any(flipsides):
                reflected_rays[flipsides,:] = this_refractedray[~reflect_here,:]
                surface_normals[flipsides,:] = -surface_normals[flipsides,:]
                n_temp = n1[flipsides]
                n1[flipsides] = n2[flipsides]
                n2[flipsides] = n_temp
                insurface[flipsides] = ~insurface[flipsides]
                still_crossing[flipsides] = np.sum(reflected_rays[flipsides,0:3]*surface_normals[flipsides,:], axis=1) <= 0
            
            # % then reflected rays
            if np.any(samesides):
                reflection_roll = np.random.randn(np.sum(reflect_here), 1)
                facet_reflection = reflection_roll.flatten() < reflectionprobs[samesides, 0]
                smooth_reflection = np.logical_and(~facet_reflection, (reflection_roll.flatten() < reflectionprobs[samesides, 1]))
                back_reflection = np.logical_and((~np.logical_or(facet_reflection, smooth_reflection)),(reflection_roll.flatten() < reflectionprobs[samesides, 2]))
                diffuse_reflection = ~np.logical_or.reduce((facet_reflection, smooth_reflection, back_reflection))

                if np.any(facet_reflection):
                    facet_ref = np.copy(samesides)
                    facet_ref[samesides] = facet_reflection
                    facet_ref_short = np.copy(reflect_here)
                    facet_ref_short[reflect_here] = facet_reflection
                    reflected_rays[facet_ref,:] = this_reflectedray[facet_ref_short,:]
                    still_crossing[facet_ref] = np.sum(reflected_rays[facet_ref,0:3]*surface_normals[facet_ref,:], axis=1) <= 0
    
                
                if np.any(smooth_reflection):
                    smooth_ref = np.copy(samesides)
                    smooth_ref[samesides] = smooth_reflection
                    [throwaway_var, theserays] = RefractionReflectionAtInterface.RefractionReflectionAtInterface(reflected_rays[smooth_ref, :], surface_normals[smooth_ref, :], n1[smooth_ref], n2[smooth_ref])
                    print("theserays smooth: " + str(np.any(np.isnan(theserays))))
                    theserays[:, 6:10] = theserays[:, 6:10] * np.matlib.repmat((reflected_rays[smooth_ref, 6] / theserays[:, 6])[:,np.newaxis], 1, 4)
                    print("theserays smooth: " + str(np.any(np.isnan(theserays))))
                    if np.any(np.isnan(theserays[:])):
                        print('whoops smooth!')
                    reflected_rays[smooth_ref, :] = theserays
                    still_crossing[smooth_ref] = False
                
                if np.any(back_reflection):
                    back_ref = np.copy(samesides)
                    back_ref[samesides] = back_reflection
                    [throwaway_var, theserays] = RefractionReflectionAtInterface.RefractionReflectionAtInterface(reflected_rays[back_ref, :], -reflected_rays[back_ref, 0:3], n1[back_ref], n2[back_ref])
                    print("theserays refl: " + str(np.any(np.isnan(theserays))))
                    theserays[:, 6:10] = theserays[:, 6:10] * np.matlib.repmat(reflected_rays[back_ref, 6] / theserays[:, 6], 1, 4)
                    print("theserays refl: " + str(np.any(np.isnan(theserays))))
                    if np.any(np.isnan(theserays[:])):
                        print('whoops back!')
                    reflected_rays[back_ref, :] = theserays
                    still_crossing[back_ref] = False
                
                if np.any(diffuse_reflection):
                    diffuse_ref = np.copy(samesides)
                    diffuse_ref[samesides] = diffuse_reflection
                    diffuse_normal = np.array(GetLambertianNormal(reflected_rays[diffuse_ref, 0:3], surface_normals[diffuse_ref,:], interface_normals[diffuse_ref,:], interface_yaxis[diffuse_ref,:])[0])
                    [throwaway_var, theserays] = RefractionReflectionAtInterface.RefractionReflectionAtInterface(reflected_rays[diffuse_ref, :], diffuse_normal, n1[diffuse_ref], n2[diffuse_ref])
                    print("theserays diff: " + str(np.any(np.isnan(theserays))))
                    theserays[:, 6:10] = theserays[:, 6:10] * np.matlib.repmat((reflected_rays[diffuse_ref, 6] / theserays[:, 6])[:,np.newaxis], 1, 4)
                    print("theserays diff: " + str(np.any(np.isnan(theserays))))
                    if np.any(np.isnan(theserays[:])):
                        print('whoops diffuse!')
                    reflected_rays[diffuse_ref, :] = theserays
                    still_crossing[diffuse_ref] = False

        still_scattering = np.logical_and(still_scattering, insurface.flatten())

        # reflected Stokes parameters nan before here
        #print("reflected check bef: " + str(reflected_rays))
        [throwaway_var, reflected_rays[still_scattering, 0:3]] = GetLambertianNormal(reflected_rays[still_scattering, 0:3], -surface_normals[still_scattering,:], interface_normals[still_scattering,:], interface_yaxis[still_scattering,:])
        #print("reflected check aft: " + str(reflected_rays))
        reflected_rays[still_scattering, 7:10] = 0
        reflected_rays[still_scattering, 6] = reflected_rays[still_scattering, 6] * reflector_params[still_scattering, 1]
        reflected_rays[still_scattering, 3:6] = np.cross(np.matlib.repmat([1, 0, 0], np.sum(still_scattering), 1), reflected_rays[still_scattering, 0:3])
        bad_polref = np.logical_and(still_scattering, (np.sum(reflected_rays[:, 3:6]**2, axis=1)<=0))
        reflected_rays[bad_polref, 3:6] = np.cross(np.matlib.repmat([0, 1, 0], np.sum(bad_polref), 1), reflected_rays[bad_polref, 0:3])
        reflected_rays[still_scattering, 3:6] = reflected_rays[still_scattering, 3:6] / np.matlib.repmat(np.abs(np.sqrt(np.sum(reflected_rays[still_scattering, 3:6]**2,1)))[:, np.newaxis],1,3)

    if np.any(np.isnan(reflected_rays[:])):
        print('whoops!')

    # %% all done!
    return reflected_rays

