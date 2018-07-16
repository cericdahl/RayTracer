% function reflected_rays = ...
%     UnifiedReflectorModel(incoming_rays, surface_normals, n1, n2, reflector_params)
%
% inputs:
%           incoming_rays   -   N-by-10 vector, where N is the number of
%                                 rays to refract.  The first three
%                                 elements are the x, y, and z components
%                                 of the forward direction of the ray (will
%                                 be normalized if not already), the 4-6
%                                 elements define the polarization
%                                 reference frame (will be projected into
%                                 plane perpendicular to ray and normalized
%                                 if not already), and the 7-10 elements
%                                 are the stokes parameters for the ray (7
%                                 gives intensity, 8-10 give polarization)
%           surface_normals -   N-by-3 vector giving the outward-pointing
%                                 normals on the incident surface
%           n1              -   1 or N element vector giving the index of
%                                 refraction of the medium the incoming ray
%                                 travels through
%           n2              -   1 or N element vector giving the index of
%                                 refraction of the medium the on the other
%                                 side of the interface
%           reflector_params-   1-by-6 or N-by-6 array giving the UNIFIED
%                                 parameters for this reflecting surface.
%                                 Columns 1:6 are:
%                                   (:,1) - surface roughness sigma,
%                                           in radians
%                                   (:,2) - reflection probability for
%                                           "transmitted" rays (always
%                                           diffuse, but subject to surface
%                                           roughness on way out)
%                                   (:,3) - lobed specular reflection
%                                           probability for "reflected"
%                                           rays
%                                   (:,4) - smooth specular reflection
%                                           probability for "reflected"
%                                           rays
%                                   (:,5) - retro reflection
%                                           probability for "reflected"
%                                           rays
% outputs:
%           reflected_rays  -   N-by-10, same format as incoming_rays, will
%                                 be normalized
%
% 8/18/16, CED

function reflected_rays = ...
    UnifiedReflectorModel(incoming_rays, surface_normals, n1, n2, reflector_params)

reflected_rays = [];

%% check inputs
if nargin<5 || size(reflector_params,2)~=5 || size(surface_normals,2)~=3 || size(incoming_rays,2)~=10 || ...
        (size(incoming_rays,1)~=size(reflector_params,1) && size(reflector_params,1)~=1) || ...
        size(incoming_rays,1)~=size(surface_normals,1) || ...
        (numel(n1)~=1 && numel(n1)~=size(incoming_rays,1)) || (numel(n2)~=1 && numel(n2)~=size(incoming_rays,1))
    disp('impropper input to UnifiedReflectorModel');
    return
end

if numel(n1)==1
    n1 = repmat(n1,size(incoming_rays,1),1);
else
    n1 = n1(:);
end

if numel(n2)==1
    n2 = repmat(n2,size(incoming_rays,1),1);
else
    n2 = n2(:);
end

if size(reflector_params,1)==1
    reflector_params = repmat(reflector_params,size(incoming_rays,1),1);
end

reflectionprobs = cumsum(reflector_params(:, 3:5), 2);
if any(reflectionprobs(:)>1)
    disp('impropper reflection parameters to UnifiedReflectorModel');
    return
end

%% normalize inputs
goodray_cut = sum(incoming_rays(:,1:3).^2,2) > 0;
if any(goodray_cut)
    incoming_rays(goodray_cut,1:3) = incoming_rays(goodray_cut,1:3) ./ ...
        repmat(abs(sqrt(sum(incoming_rays(goodray_cut,1:3).^2,2))),1,3);
end

goodsurface_cut = sum(surface_normals.^2,2) > 0;
if any(goodsurface_cut)
    surface_normals(goodsurface_cut,:) = surface_normals(goodsurface_cut,:) ./ ...
        repmat(abs(sqrt(sum(surface_normals(goodsurface_cut,:).^2,2))),1,3);
end

incoming_rays(:,4:6) = incoming_rays(:,4:6) - ...
    repmat(sum(incoming_rays(:,4:6).*incoming_rays(:,1:3),2),1,3) .* incoming_rays(:,1:3);
goodpolarization_cut = sum(incoming_rays(:,4:6).^2,2) > 0;
if any(goodpolarization_cut)
    incoming_rays(goodpolarization_cut,4:6) = incoming_rays(goodpolarization_cut,4:6) ./ ...
        repmat(abs(sqrt(sum(incoming_rays(goodpolarization_cut,4:6).^2,2))),1,3);
end

%% set defaults
reflected_rays = incoming_rays;

%% find interface normals
cos_incident_angle = sum(-incoming_rays(:,1:3).*surface_normals, 2);
goodhit_cut = cos_incident_angle > 0;

interface_normals = cross(-incoming_rays(:,1:3), surface_normals, 2);
sin_incident_angle = abs(sqrt(sum(interface_normals.^2, 2)));
goodinterface_cut = sin_incident_angle > 0;
if any(goodinterface_cut)
    interface_normals(goodinterface_cut,:) = interface_normals(goodinterface_cut,:) ./ ...
        repmat(sin_incident_angle(goodinterface_cut),1,3);
end
if any(~goodinterface_cut)
    tmp_inorms = cross(repmat([1,0,0], sum(~goodinterface_cut), 1), ...
        surface_normals(~goodinterface_cut, :), 2);
    tmp_norms_good = sum(tmp_inorms.^2, 2)>0;
    tmp_inorms2 = cross(repmat([0,1,0], sum(~goodinterface_cut), 1), ...
        surface_normals(~goodinterface_cut, :), 2);
    tmp_inorms(~tmp_norms_good,:) = tmp_inorms2(~tmp_norms_good,:);
    interface_normals(~goodinterface_cut,:) = tmp_inorms ./ ...
        repmat(abs(sqrt(sum(tmp_inorms.^2, 2))), 1, 3);
end

%% and complete the local bases of (interface_normal, interface_yaxis, surface_normal)
interface_yaxis = cross(surface_normals, interface_normals, 2);

%% Now the fun starts
still_scattering = goodhit_cut;
insurface = false(size(incoming_rays,1),1);

% loop until all the rays have been reflected
while any(still_scattering)
    
    % loop until all the rays have finished their interaction with the
    % rough dielectric boundary
    still_crossing = still_scattering;
    while any(still_crossing)
        % First find a microfacet
        facet_normals = GetFacetNormal(reflected_rays(still_crossing, 1:3), ...
            surface_normals(still_crossing,:), interface_normals(still_crossing,:), ...
            interface_yaxis(still_crossing,:), reflector_params(still_crossing,1));

        % then get refrations/reflections off of this facet
        [this_refractedray, this_reflectedray] = RefractionReflectionAtInterface(reflected_rays(still_crossing, :), ...
            facet_normals, n1(still_crossing), n2(still_crossing));
    
        % roll dice to see whether we're following the refraction or the
        % reflection
        reflect_here = rand(size(facet_normals,1),1) > ...
            (this_refractedray(:,7) ./ (this_refractedray(:,7) + this_reflectedray(:,7)));
        
        if any(isnan(this_refractedray(:)))
            disp('whoops!');
        end
        if any(isnan(this_reflectedray(:)))
            disp('whoops!');
        end
        
        % renormalize refracted/reflected intensities to the initial ray
        this_refractedray(:,7:10) = this_refractedray(:,7:10) .* ...
            repmat(reflected_rays(still_crossing,7)./this_refractedray(:,7), 1, 4);
        this_reflectedray(:,7:10) = this_reflectedray(:,7:10) .* ...
            repmat(reflected_rays(still_crossing,7)./this_reflectedray(:,7), 1, 4);
        
        if any(isnan(this_refractedray(~reflect_here, 7)))
            disp('whoops!');
        end
        if any(isnan(this_reflectedray(reflect_here, 7)))
            disp('whoops!');
        end
        
        flipsides = still_crossing;
        flipsides(still_crossing) = ~reflect_here;
        samesides = still_crossing;
        samesides(still_crossing) = reflect_here;

        % first handle refracted rays
        if any(flipsides)
            reflected_rays(flipsides,:) = this_refractedray(~reflect_here,:);
            surface_normals(flipsides,:) = -surface_normals(flipsides,:);
            n_temp = n1(flipsides);
            n1(flipsides) = n2(flipsides);
            n2(flipsides) = n_temp;
            insurface(flipsides) = ~insurface(flipsides);
            still_crossing(flipsides) = sum(reflected_rays(flipsides,1:3).*surface_normals(flipsides,:), 2) <= 0;
        end
        
        % then reflected rays
        if any(samesides)
            reflection_roll = rand(sum(reflect_here),1);
            facet_reflection = reflection_roll < reflectionprobs(samesides, 1);
            smooth_reflection = ~facet_reflection & (reflection_roll < reflectionprobs(samesides, 2));
            back_reflection = ~(facet_reflection | smooth_reflection) & (reflection_roll < reflectionprobs(samesides, 3));
            diffuse_reflection = ~(facet_reflection  | smooth_reflection | back_reflection);

            if any(facet_reflection)
                facet_ref = samesides;
                facet_ref(samesides) = facet_reflection;
                facet_ref_short = reflect_here;
                facet_ref_short(reflect_here) = facet_reflection;
                reflected_rays(facet_ref,:) = this_reflectedray(facet_ref_short,:);
                still_crossing(facet_ref) = sum(reflected_rays(facet_ref,1:3).*surface_normals(facet_ref,:), 2) <= 0;
            end
            
            if any(smooth_reflection)
                smooth_ref = samesides;
                smooth_ref(samesides) = smooth_reflection;
                [~, theserays] = RefractionReflectionAtInterface(reflected_rays(smooth_ref, :), ...
                    surface_normals(smooth_ref, :), n1(smooth_ref), n2(smooth_ref));
                theserays(:, 7:10) = theserays(:, 7:10) .* ...
                    repmat(reflected_rays(smooth_ref, 7) ./ theserays(:, 7), 1, 4);
                if any(isnan(theserays(:)))
                    disp('whoops!');
                end
                reflected_rays(smooth_ref, :) = theserays;
                still_crossing(smooth_ref) = false;
            end
            
            if any(back_reflection)
                back_ref = samesides;
                back_ref(samesides) = back_reflection;
                [~, theserays] = RefractionReflectionAtInterface(reflected_rays(back_ref, :), ...
                    -reflected_rays(back_ref, 1:3), n1(back_ref), n2(back_ref));
                theserays(:, 7:10) = theserays(:, 7:10) .* ...
                    repmat(reflected_rays(back_ref, 7) ./ theserays(:, 7), 1, 4);
                if any(isnan(theserays(:)))
                    disp('whoops!');
                end
                reflected_rays(back_ref, :) = theserays;
                still_crossing(back_ref) = false;
            end
            
            if any(diffuse_reflection)
                diffuse_ref = samesides;
                diffuse_ref(samesides) = diffuse_reflection;
                diffuse_normal = GetLambertianNormal(reflected_rays(diffuse_ref, 1:3), surface_normals(diffuse_ref,:), ...
                    interface_normals(diffuse_ref,:), interface_yaxis(diffuse_ref,:));
                [~, theserays] = RefractionReflectionAtInterface(reflected_rays(diffuse_ref, :), ...
                    diffuse_normal, n1(diffuse_ref), n2(diffuse_ref));
                theserays(:, 7:10) = theserays(:, 7:10) .* ...
                    repmat(reflected_rays(diffuse_ref, 7) ./ theserays(:, 7), 1, 4);
                if any(isnan(theserays(:)))
                    disp('whoops!');
                end
                reflected_rays(diffuse_ref, :) = theserays;
                still_crossing(diffuse_ref) = false;
            end
        end
    end
    
    still_scattering = still_scattering & insurface;
    
    [~, reflected_rays(still_scattering, 1:3)] = GetLambertianNormal(reflected_rays(still_scattering, 1:3), ...
        -surface_normals(still_scattering,:), interface_normals(still_scattering,:), interface_yaxis(still_scattering,:));
    reflected_rays(still_scattering, 8:10) = 0;
    reflected_rays(still_scattering, 7) = reflected_rays(still_scattering, 7) .* reflector_params(still_scattering, 2);
    reflected_rays(still_scattering, 4:6) = cross(repmat([1, 0, 0], sum(still_scattering), 1), reflected_rays(still_scattering, 1:3));
    bad_polref = still_scattering & sum(reflected_rays(:, 4:6).^2, 2)<=0;
    reflected_rays(bad_polref, 4:6) = cross(repmat([0, 1, 0], sum(bad_polref), 1), reflected_rays(bad_polref, 1:3));
    reflected_rays(still_scattering, 4:6) = reflected_rays(still_scattering, 4:6) ./ ...
        repmat(abs(sqrt(sum(reflected_rays(still_scattering, 4:6).^2,2))),1,3);
    
end

if any(isnan(reflected_rays(:)))
    disp('whoops!');
end

%% all done!
return

function facet_normal = GetFacetNormal(indir, s_norm, s_x, s_y, sig_a)
    % Same implementation here as in Geant4
    facet_normal = s_norm;
    facets_set = sig_a==0;
    while ~all(facets_set)
        these_sig_a = sig_a(~facets_set);
        thetas = abs(these_sig_a .* randn(sum(~facets_set), 1));
        
        costhetas = cos(thetas);
        sinthetas = sin(thetas);
        
        out_of_range = thetas >= (.5*pi);
        % next line changes prob distribution from gauss to sin*gauss --
        % the these_sig_a*4 makes this approximate out to 4-sigma, same as
        % in G4OpBoundaryProcess::GetFacetNormal
        fail_prob_jacob = rand(size(thetas)).*min(these_sig_a*4,1) > sinthetas;
        
        phis = 2*pi*rand(size(thetas));
        
        facet_normal(~facets_set, :) = s_norm(~facets_set,:) .* repmat(costhetas, 1, 3) + ...
            s_x(~facets_set,:) .* repmat(sinthetas, 1, 3) .* repmat(cos(phis), 1, 3) + ...
            s_y(~facets_set,:) .* repmat(sinthetas, 1, 3) .* repmat(sin(phis), 1, 3);
        
        wrongside = sum(indir(~facets_set,:).*facet_normal(~facets_set,:), 2)>=0;
        
        facets_set(~facets_set) = ~(out_of_range | fail_prob_jacob | wrongside);
    end
    
return


function [facet_normal, outdir] = GetLambertianNormal(indir, s_norm, s_x, s_y)

    out_costheta = sqrt(rand(size(indir,1), 1)); % diffuse ~= isotropic
    out_sintheta = sqrt(1 - out_costheta.^2);
    out_phi = 2*pi*rand(size(out_costheta));
    outdir = s_norm .* repmat(out_costheta, 1, 3) + ...
        s_x .* repmat(out_sintheta, 1, 3) .* repmat(cos(out_phi), 1, 3) + ...
        s_y .* repmat(out_sintheta, 1, 3) .* repmat(sin(out_phi), 1, 3);
    facet_normal = outdir - indir;
    facet_normal = facet_normal ./ repmat(abs(sqrt(sum(facet_normal.^2,2))),1,3);
    
return