% function [refracted_rays reflected_rays] = ...
%     RefractionReflectionAtInterface(incoming_rays, surface_normals, n1, n2, tir_handling)
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
%           tir_handling    -   This determines what the refracted_rays
%                                 output is in the case of total internal
%                                 reflection.  The default (-1) gives a
%                                 refracted ray tangent to the surface with
%                                 zero intensity.  Any value >=0 will give
%                                 a ray with the same direction and
%                                 polarization as the reflected ray, with
%                                 intensity equal to the reflected
%                                 intensity times tir_handling.  This lets
%                                 you treat tir-rays like refracted rays,
%                                 which can be handy in geometry sims.
% outputs:
%           refracted_rays  -   N-by-10, same format as incoming_rays, will
%                                 be normalized
%           reflected_rays  -   N-by-10, same format as incoming_rays, will
%                                 be normalized
%
% 12/15/09, CED

function [refracted_rays, reflected_rays] = ...
    RefractionReflectionAtInterface(incoming_rays, surface_normals, n1, n2, tir_handling)

refracted_rays = [];
reflected_rays = [];

%% check inputs
if nargin<5 || isempty(tir_handling)
    tir_handling = -1;
end

if nargin<4 || size(surface_normals,2)~=3 || size(incoming_rays,2)~=10 || size(incoming_rays,1)~=size(surface_normals,1) || ...
        (numel(n1)~=1 && numel(n1)~=size(incoming_rays,1)) || (numel(n2)~=1 && numel(n2)~=size(incoming_rays,1)) || ...
        (numel(tir_handling)~=1 && numel(tir_handling)~=size(incoming_rays,1))
    disp('impropper input to RefractionReflactionAtInterface');
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

if numel(tir_handling)==1
    tir_handling = repmat(tir_handling,size(incoming_rays,1),1);
else
    tir_handling = tir_handling(:);
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
refracted_rays = incoming_rays;
reflected_rays = incoming_rays;
reflected_rays(:,1:3) = -reflected_rays(:,1:3);
reflected_rays(:,7:10) = 0;

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

%% rotate stokes parameters to new basis
c_rot = sum(interface_normals .* incoming_rays(:,4:6), 2);
s_rot = sum(cross(interface_normals,incoming_rays(:,4:6),2) .* ...
    incoming_rays(:,1:3), 2);
c2_rot = c_rot.^2 - s_rot.^2;
s2_rot = 2.*c_rot.*s_rot;
% for scalar c2_rot and s2_rot, we would write:
%     rotmat = [c2_rot s2_rot ; -s2_rot c2_rot];
% but the transformation is different for every ray, so.......
if any(goodinterface_cut)
    incoming_rays(goodinterface_cut,4:6) = interface_normals(goodinterface_cut,:);
% again, if we had a 2x2 rotmat, we would write:
%     incoming_rays(goodinterface_cut,8:9) = incoming_rays(goodinterface_cut,8:9) * rotmat;
% but instead we do it one element at a time
    old_polarization = incoming_rays(goodinterface_cut,8:9);
    incoming_rays(goodinterface_cut,8) = old_polarization(:,1).*c2_rot(goodinterface_cut) - old_polarization(:,2).*s2_rot(goodinterface_cut);
    incoming_rays(goodinterface_cut,9) = old_polarization(:,1).*s2_rot(goodinterface_cut) + old_polarization(:,2).*c2_rot(goodinterface_cut);
    refracted_rays(goodinterface_cut,4:10) = incoming_rays(goodinterface_cut,4:10);
    reflected_rays(goodinterface_cut,4:6) = incoming_rays(goodinterface_cut,4:6);
end

%% find complex amplitudes for incoming rays in new basis
% We break the incoming ray into 3 polarized planes waves -- one in the
% polarization direction indicated by the stokes parameters, one in the
% interface_normal, and one in the interface-plane -- the latter two make
% up the unpolarized portion of the incoming ray
amplitudes = zeros(size(incoming_rays,1),3,2);
polarized_intensity = abs(sqrt(sum(incoming_rays(:,8:10).^2,2)));
amplitudes(:,1,1) = abs(sqrt(.5 .* (polarized_intensity + incoming_rays(:,8))));
amplitudes(:,1,2) = abs(sqrt(.5 .* (polarized_intensity - incoming_rays(:,8)))) .* ...
    exp(1i.*atan2(incoming_rays(:,10),incoming_rays(:,9)));
amplitudes(:,2,1) = abs(sqrt(.5 .* (incoming_rays(:,7) - polarized_intensity)));
amplitudes(:,3,2) = amplitudes(:,2,1);

%% calculated reflected and refracted amplitudes
sin_refracted_angle = sin_incident_angle .* n1 ./ n2;
cos_refracted_angle = sqrt(1 - sin_refracted_angle.^2);

rs = (n1.*cos_incident_angle - n2.*cos_refracted_angle) ./ ...
    (n1.*cos_incident_angle + n2.*cos_refracted_angle);
rp = -(n1.*cos_refracted_angle - n2.*cos_incident_angle) ./ ...
    (n1.*cos_refracted_angle + n2.*cos_incident_angle);

rs(n2==inf | n2==-inf) = -1;
rp(n2==inf | n2==-inf) = 1;

ts = abs(sqrt(1-conj(rs).*rs));
tp = abs(sqrt(1-conj(rp).*rp));
refracted_amplitudes = amplitudes .* repmat(reshape([ts tp],[],1,2),[1 3 1]);
reflected_amplitudes = amplitudes .* repmat(reshape([rs rp],[],1,2),[1 3 1]);

%% get back to stokes parameters
if any(goodhit_cut)
    refracted_rays(goodhit_cut,7) = sum(sum(conj(refracted_amplitudes(goodhit_cut,:,:)) .* ...
        refracted_amplitudes(goodhit_cut,:,:),3),2);
    refracted_rays(goodhit_cut,8) = -sum(diff(conj(refracted_amplitudes(goodhit_cut,:,:)) .* ...
        refracted_amplitudes(goodhit_cut,:,:),1,3),2);
    refracted_rays(goodhit_cut,9) = sum(real(2 .* conj(refracted_amplitudes(goodhit_cut,:,1)) .* ...
        refracted_amplitudes(goodhit_cut,:,2)),2);
    refracted_rays(goodhit_cut,10) = sum(imag(2 .* conj(refracted_amplitudes(goodhit_cut,:,1)) .* ...
        refracted_amplitudes(goodhit_cut,:,2)),2);

    reflected_rays(goodhit_cut,7) = sum(sum(conj(reflected_amplitudes(goodhit_cut,:,:)) .* ...
        reflected_amplitudes(goodhit_cut,:,:),3),2);
    reflected_rays(goodhit_cut,8) = -sum(diff(conj(reflected_amplitudes(goodhit_cut,:,:)) .* ...
        reflected_amplitudes(goodhit_cut,:,:),1,3),2);
    reflected_rays(goodhit_cut,9) = sum(real(2 .* conj(reflected_amplitudes(goodhit_cut,:,1)) .* ...
        reflected_amplitudes(goodhit_cut,:,2)),2);
    reflected_rays(goodhit_cut,10) = sum(imag(2 .* conj(reflected_amplitudes(goodhit_cut,:,1)) .* ...
        reflected_amplitudes(goodhit_cut,:,2)),2);
end

%% consider surface_normal to be -x, interface_normal to be +z
new_yaxis = cross(surface_normals, interface_normals, 2);
goodcut = goodhit_cut & (sum(new_yaxis.^2,2) > 0);

if any(goodcut)
    new_yaxis(goodcut,:) = new_yaxis(goodcut,:) ./ ...
        repmat(abs(sqrt(sum(new_yaxis(goodcut,:).^2,2))),1,3);
end

%% calculate reflected ray direction for non-normal incidence
if any(goodcut)
    reflected_rays(goodcut,1:3) = ...
        repmat(cos_incident_angle(goodcut),1,3) .* surface_normals(goodcut,:) - ...
        repmat(sin_incident_angle(goodcut),1,3) .* new_yaxis(goodcut,:);
end

%% calculated refracted ray direction
total_internal_reflection_cut = goodcut & sin_refracted_angle >= 1;
refracted_cut = goodcut & ~total_internal_reflection_cut;

if any(refracted_cut)
    refracted_rays(refracted_cut, 1:3) = ...
        -repmat(cos_refracted_angle(refracted_cut),1,3) .* surface_normals(refracted_cut,:) - ...
        repmat(sin_refracted_angle(refracted_cut),1,3) .* new_yaxis(refracted_cut,:);
end

if any(total_internal_reflection_cut)
    if any(total_internal_reflection_cut & (tir_handling<0))
        refracted_rays(total_internal_reflection_cut & (tir_handling<0), 1:3) = -new_yaxis(total_internal_reflection_cut & (tir_handling<0),:);
    end
    if any(total_internal_reflection_cut & (tir_handling>=0))
        refracted_rays(total_internal_reflection_cut & (tir_handling>=0), 1:6) = reflected_rays(total_internal_reflection_cut & (tir_handling>=0),1:6);
        refracted_rays(total_internal_reflection_cut & (tir_handling>=0), 7:10) = reflected_rays(total_internal_reflection_cut & (tir_handling>=0),7:10) .* ...
            repmat(tir_handling(total_internal_reflection_cut & (tir_handling>=0)),1,4);
    end
end

%% all done!
return
