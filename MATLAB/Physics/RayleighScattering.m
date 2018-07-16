% function scattered_rays = ...
%     RayleighScattering(incoming_rays)
%
% This function Rayleigh scatters each of the incoming rays.  The output
% direction is randomized according to the Rayleigh scattering distribution
% (dependent on the polarization of each incoming ray) and the output
% polarization (i.e. Stokes parameters) is fixed based on the incoming
% polarization and scattering direction.  Intensity of the ray is
% unchanged.  For a nice description of Rayleigh scattering in terms of
% Stokes parameters, see https://www.astro.umd.edu/~jph/notes3.pdf
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
%           azimuth_precision - The plane of scattering is determined by
%                                 a dice roll to a lookup table of length
%                                 (azimuth_precision) with equally-spaced
%                                 entries in cumulative probability (and
%                                 linear interpolation between entries).
%                                 Default value is 1e5.
% outputs:
%           scattered_rays  -   N-by-10, same format as incoming_rays, will
%                                 be normalized
%
%
% 8/18/16, CED

function scattered_rays = RayleighScattering(incoming_rays, azimuth_precision)

persistent rayleigh_azimuth

%% set defaults
scattered_rays = [];

if nargin<2 || isempty(azimuth_precision)
    azimuth_precision = 1e5;
end

%% check inputs
if nargin<1 || size(incoming_rays,2)~=10
    disp('impropper input to RayleighScattering');
    return
end

%% initialize rayleigh_azimuth
if isempty(rayleigh_azimuth) || length(rayleigh_azimuth) ~= azimuth_precision
    cumd = linspace(0, 1, azimuth_precision+1);
    phid = linspace(0, 2*pi, azimuth_precision*10);
    cumd_phi = (phid - .25*sin(2*phid)) / (2*pi);
    rayleigh_azimuth = interp1(cumd_phi, phid, cumd)';
end

%% normalize inputs
goodray_cut = sum(incoming_rays(:,1:3).^2,2) > 0;
if any(goodray_cut)
    incoming_rays(goodray_cut,1:3) = incoming_rays(goodray_cut,1:3) ./ ...
        repmat(abs(sqrt(sum(incoming_rays(goodray_cut,1:3).^2,2))),1,3);
end

incoming_rays(:,4:6) = incoming_rays(:,4:6) - ...
    repmat(sum(incoming_rays(:,4:6).*incoming_rays(:,1:3),2),1,3) .* incoming_rays(:,1:3);
goodpolarization_cut = sum(incoming_rays(:,4:6).^2,2) > 0;
if any(goodpolarization_cut)
    incoming_rays(goodpolarization_cut,4:6) = incoming_rays(goodpolarization_cut,4:6) ./ ...
        repmat(abs(sqrt(sum(incoming_rays(goodpolarization_cut,4:6).^2,2))),1,3);
end

%% roll dice to determine scattering angle
scatter_dice = rand(size(incoming_rays,1), 3);
% column 1 for the lin-polarized vs not-lin-polarized bit
% column 2 for phi
% column 3 for theta

%% Determine scattering plane relative to polarization reference
linpol_intensity = abs(sqrt(sum(incoming_rays(:,8:9).^2,2)));
linpol_frac = linpol_intensity ./ incoming_rays(:,7);
linpol_scatter = linpol_frac > scatter_dice(:, 1);

scatter_phi = 2*pi*scatter_dice(:, 2);

if any(linpol_scatter)
    linpol_angle = .5 * atan2(incoming_rays(linpol_scatter,9), incoming_rays(linpol_scatter,8));
    ix = floor(azimuth_precision*scatter_dice(linpol_scatter, 2)) + 1;
    ix = min(ix, azimuth_precision);
    ix_frac = azimuth_precision*scatter_dice(linpol_scatter, 2) - (ix-1);
    scatter_phi(linpol_scatter) = linpol_angle + ...
        rayleigh_azimuth(ix) + ...
        ix_frac.*(rayleigh_azimuth(ix+1) - rayleigh_azimuth(ix));
end

%% Now rotate stokes parameters to the scattering plane
c_rot = cos(scatter_phi);
s_rot = sin(scatter_phi);
c2_rot = c_rot.^2 - s_rot.^2;
s2_rot = 2.*c_rot.*s_rot;

old_polarization = incoming_rays(:,8:9);
incoming_rays(:,8) = old_polarization(:,1).*c2_rot + old_polarization(:,2).*s2_rot;
incoming_rays(:,9) = -old_polarization(:,1).*s2_rot + old_polarization(:,2).*c2_rot;

old_reference = incoming_rays(:,4:6);
old_reference_perp = cross(incoming_rays(:, 1:3), incoming_rays(:, 4:6));
incoming_rays(:, 4:6) = old_reference .* repmat(c_rot, 1, 3) + ...
    old_reference_perp .* repmat(s_rot, 1, 3);

new_reference_perp = cross(incoming_rays(:, 1:3), incoming_rays(:, 4:6));

%% Now determine scattering angle
a = incoming_rays(:, 8) ./ incoming_rays(:, 7); % -1 < a < 1
% now cos_theta is the solution of a cubic polynomial with coefficients
% determined by a and scatter_dice(:, 3)
z = (2-a) .* (1 - 2*scatter_dice(:, 3));
s = sqrt(z.^2 + (1-a).^3./(1+a));
A = sign(z+s) .* abs(z + s).^(1/3);
B = sign(z-s) .* abs(z - s).^(1/3);
cos_theta = (A + B) .* (1 + a).^(-1/3);
sin_theta = sqrt(1 - cos_theta.^2);

%% now determine new direction
scattered_rays = incoming_rays;
scattered_rays(:,1:3) = incoming_rays(:,1:3) .* repmat(cos_theta, 1, 3) + ...
    old_reference .* repmat(sin_theta, 1, 3) .* repmat(c_rot, 1, 3) + ...
    old_reference_perp .* repmat(sin_theta, 1, 3) .* repmat(s_rot, 1, 3);
scattered_rays(:, 4:6) = cross(new_reference_perp, scattered_rays(:, 1:3));

%% and new ray polarization
c2theta = cos_theta.^2;

R11 = c2theta + 1;
R12 = c2theta - 1;

rnorm = incoming_rays(:, 7) ./ (incoming_rays(:, 7) .* R11 + incoming_rays(:, 8) .* R12);
scattered_rays(:, 8) = (incoming_rays(:, 7) .* R12 + incoming_rays(:, 8) .* R11) .* rnorm;

R33 = 2*cos_theta .* rnorm;

scattered_rays(:, 9) = incoming_rays(:, 9) .* R33;
scattered_rays(:, 10) = incoming_rays(:, 10) .* R33;

