% function surface_list = Create2LGeometry(geospecs)
%
% This function creates a structure array of surfaces to be used by
% RayTracer.  Follow this architecture to create any geometry you like.
%
% Each surface has six fields.  The first (intersect_function) defines the
% geometry, and is a function handle to an anonymous function, calling a
% RayToXXXXX function with the appropriate geometry inputs.  For example, 
%
%  @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], 10) 
%
% defines a cylinder on the z-axis with radius 10.  See all RayToXXXXX
% functions in the RayTracing directory for other possible shapes (and
% create your own if you desire).
%
% The second field (inbounds_function) defines the bounds of the surface,
% and is a function handle to an anonymous function that inputs an N-by-3-by-M 
% array and outputs an N-by-M logical.  It can be assumed that all input
% points are on the surface defined by intersect_function, giving true if
% the input point is contained within the bounds, false otherwise.  For
% example,
%
%  @(p)(reshape( ...
%      (p(:,3,:)>20) & (p(:,3,:)<80) & (atan2(p(:,2,:),p(:,1,:))>0), ...
%      size(p,1), [] ));
%
% would cut the above cylinder in half along the xz plane and truncate it 
% at 20 and 80 in the z-coordinate.  (Generically, p is an N-by-3-by-M
% matrix, and the output of the function should be an N-by-M boolean.)
%
% The third and fourth fields are n_outside and n_inside, and give the
% indices of refraction on the two sides of the surface.  What is 'inside'
% and 'outside' is defined in the RayToXXXXX function -- for spheres and
% cylinders, this is obvious, for planes less so.  See the documentation
% for each RayToXXXXX function for details.  Also, setting n to inf makes
% that side a perfect conductor (in terms of calculating reflection and
% polarization, at least).
%
% The fifth field is surface type, and may be 'normal', 'diffuse', or
% 'retro'.  For 'diffuse' and 'retro' surfaces, the normal direction
% returned by the intersect_function is replaced by a random direction
% within pi/4 of normal or the reverse of the incoming ray, respectively.
%
% The sixth field is an absorption coefficient -- RayTracer will multiply
% the intensity of both the reflected and refracted rays coming from this
% surface by 1-absorption.
%
% 12/16/09, CED

function [surface_list rays ray_startingpoints pixels test_pixels] = Create2LGeometry(geospecs)

if nargin<1 || isempty(geospecs) || ~isstruct(geospecs)
    geospecs = struct();
end

%% define surface_list structure
surface_list = struct( ...
    'description', {}, ...
    'intersect_function', {}, ...
    'inbounds_function', {}, ...
    'n_outside', {}, ...
    'n_inside', {}, ...
    'surface_type', {}, ...
    'absorption', {});

rays = cell(1,2);
ray_startingpoints = cell(1,2);
pixels = cell(1,2);
test_pixels = cell(1,2);

%% define default values for everything that may be included in geospecs
% indices of refraction
n_CF3I = 1.31;
n_H2O = 1.33;
n_quartz = 1.458;
n_glycol = 1.434;
n_air = 1.00;
n_glass = 1.52;

% focal lengths in cm
cam0_focallength = .5;
cam1_focallength = .5;

% positions in cm, x=0, z=0 at center of hemisphere, y=0 at window-air interface
cam0_x = -4;
cam0_y = -3;
cam0_z = 2.5;
cam1_x = 4;
cam1_y = -3;
cam1_z = 2.5;

% jar dimensions in cm
jar_cylthick = .25; % thickness of cylinder wall
jar_axthick = .25; % thickness of sphere wall at apex
jar_cylrad = 7.4464; % outer radius of cylinder
jar_axrad = 8.3954; % outer radius of sphere (along cylinder axis)

% cf3i volume
cf3i_mass = 4048;
cf3i_density = 2;

% jar rotation in degrees (roll about +z, pitch about +y, yaw about +z)
jar_pitch = 5; % angle between axis and +z
jar_yaw = 90;
jar_roll = -90;

% cam rotation in degrees (roll about +y, pitch about +x, yaw about +z)
cam0_pitch = -3;
cam0_yaw = -10;
cam0_roll = 0;
cam1_pitch = -3;
cam1_yaw = 10;
cam1_roll = 0;

% cam distortion (unitless, adds quadratic term to ray direction)
cam0_distortion = 0;
cam1_distortion = 0;

% window dimensions in cm
window_inside = -(.5*11.938 + 2.8)*2.54 - .254; % -horizontal distance from hemisphere center
window_thickness = .9*2.54;

% fiducial marks in cm
fid_mark_z1 = 8.5; % distance from hemisphere apex
fid_mark_z2 = 13; % distance from hemisphere apex
fid_mark_rphi = 16; % distance around circumference to back points
fid_mark_length = .5; % length from center of mark to one crosshair
fid_mark_pen = .1; % half-width of crosshair line

% surface test marks in cm, degrees
surface_test_cyl_z = 1; % z=0 at hemisphere center
surface_test_cyl_phi = 180; % phi = 0 at -y (towards camera)
surface_test_sph_z = -1; % z=0 at hemisphere center
surface_test_sph_phi = 180; % phi = 0 at -y (towards camera)
testmark_radius = .2; % dot radius

lens_type = 'theta';

%% overwrite default values with values in geospecs
fnames = fieldnames(geospecs);
for fn=1:length(fnames)
    if ~isempty(geospecs.(fnames{fn}))
        eval([fnames{fn} ' = geospecs.(fnames{fn});']);
    end
end

%% a few fixed distances and derived quantities
cam_pixel_pitch = .00099;
cam_resolution = [491 656];

if 0
    cam_ccd_size = cam_pixel_pitch * (cam_resolution - 1);
    cam_resolution = [100 134];
    cam_pixel_pitch = cam_ccd_size ./ (cam_resolution - 1);
end

cam0_position = [cam0_x cam0_y+window_inside-window_thickness cam0_z];
cam1_position = [cam1_x cam1_y+window_inside-window_thickness cam1_z];

jar_rotate = [jar_pitch jar_yaw jar_roll]*pi/180;
cam0_rotate = [cam0_pitch cam0_yaw cam0_roll]*pi/180;
cam1_rotate = [cam1_pitch cam1_yaw cam1_roll]*pi/180;

jar_rotmat = [cos(jar_rotate(2)) -sin(jar_rotate(2)) 0 ; sin(jar_rotate(2)) cos(jar_rotate(2)) 0 ; 0 0 1] * ...
    [cos(jar_rotate(1)) 0 sin(jar_rotate(1)) ; 0 1 0 ; -sin(jar_rotate(1)) 0 cos(jar_rotate(1))] * ...
    [cos(jar_rotate(3)) -sin(jar_rotate(3)) 0 ; sin(jar_rotate(3)) cos(jar_rotate(3)) 0 ; 0 0 1];

fid_mark_phi = fid_mark_rphi / jar_cylrad;
fid_mark_lengthphi = fid_mark_length / jar_cylrad;
fid_mark_penphi = fid_mark_pen / jar_cylrad;
fid_mark_z = [fid_mark_z1 fid_mark_z2] - jar_axrad;

test_cyl = [0 0 surface_test_cyl_z] + ...
    (jar_cylrad-jar_cylthick)*[sin(surface_test_cyl_phi*pi/180) -cos(surface_test_cyl_phi*pi/180) 0];
test_cyl = test_cyl*(jar_rotmat');
test_sph = [0 0 surface_test_sph_z] + ...
    sqrt(1 - (surface_test_sph_z/(jar_axrad-jar_axthick))^2)*(jar_cylrad-jar_cylthick)*[sin(surface_test_sph_phi*pi/180) -cos(surface_test_sph_phi*pi/180) 0];
test_sph = test_sph*(jar_rotmat');

window_outside = window_inside - window_thickness; 

cf3i_volume = cf3i_mass / cf3i_density;
hemi_volume = (2/3)*pi*(jar_cylrad-jar_cylthick)^2*(jar_axrad-jar_axthick) + ...
    pi*(jar_cylrad-jar_cylthick)^3*tan(jar_rotate(1));
liquid_level_cyl = (cf3i_volume - hemi_volume) / (pi*(jar_cylrad-jar_cylthick)^2);
liquid_level = liquid_level_cyl*cos(jar_rotate(1)) + (jar_cylrad-jar_cylthick)*sin(jar_rotate(1));


%% Clausius-Mosotti
% clausius_mosotti_n = @(rho, n0, rho0)(sqrt((n0.^2 .* (rho0 + 2.*rho) + 2 .* (rho0 - rho))./(n0.^2 .* (rho0 - rho) + (2.*rho0 + rho))));

%% surface list

surface_list(end+1).description = 'inside surface of quartz cylinder below water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1]*(jar_rotmat'), jar_cylrad-jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<liquid_level) & ...
    (sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2) > 0) & ...
    (sum((p - repmat(test_cyl,[size(p,1),1,size(p,3)])).^2,2)>testmark_radius^2), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz cylinder above water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1]*(jar_rotmat'), jar_cylrad-jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=liquid_level) & ...
    (sum((p - repmat(test_cyl,[size(p,1),1,size(p,3)])).^2,2)>testmark_radius^2), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_H2O;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1]*(jar_rotmat'), jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2) > 0) & ...
    ~( ( ( (abs(reshape(abs(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) )), size(p,1), 1, []) - fid_mark_phi)<fid_mark_penphi) | ...
    (abs(reshape(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) ), size(p,1), 1, []) )<fid_mark_penphi) ) & ...
    ( (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(1))<fid_mark_length) | ...
    (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(2))<fid_mark_length) ) ) | ...
    ( ( (abs(reshape(abs(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) )), size(p,1), 1, []) - fid_mark_phi)<fid_mark_lengthphi) | ...
    (abs(reshape(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) ), size(p,1), 1, []) )<fid_mark_lengthphi) ) & ...
    ( (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(1))<fid_mark_pen) | ...
    (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(2))<fid_mark_pen) ) ) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    jar_rotmat*[(jar_cylrad-jar_cylthick)^-2 0 0 ; 0 (jar_cylrad-jar_cylthick)^-2 0 ; 0 0 (jar_axrad-jar_axthick)^-2 ]*jar_rotmat', ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2) <= 0) & ...
    (sum((p - repmat(test_sph,[size(p,1),1,size(p,3)])).^2,2)>testmark_radius^2), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    jar_rotmat*[jar_cylrad^-2 0 0 ; 0 jar_cylrad^-2 0 ; 0 0 jar_axrad^-2 ]*jar_rotmat', ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2) <= 0), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'CF3I - water interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 liquid_level], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( sum( ( ...
    p - ...
    ( repmat(sum(p .* repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]), 2), [1 3 1]) .* ...
    repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]) ) ...
    ).^2,2) < ((jar_cylrad-jar_cylthick)^2), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_H2O;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'glass - glycol interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 window_inside 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(true(size(p,1),size(p,3)));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'glass - air interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 window_outside 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(true(size(p,1),size(p,3)));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'fiducial marks';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1]*(jar_rotmat'), jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2) > 0) & ...
    ( ( ( (abs(reshape(abs(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) )), size(p,1), 1, []) - fid_mark_phi)<fid_mark_penphi) | ...
    (abs(reshape(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) ), size(p,1), 1, []) )<fid_mark_penphi) ) & ...
    ( (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(1))<fid_mark_length) | ...
    (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(2))<fid_mark_length) ) ) | ...
    ( ( (abs(reshape(abs(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) )), size(p,1), 1, []) - fid_mark_phi)<fid_mark_lengthphi) | ...
    (abs(reshape(atan2( real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([1 0 0],size(p,1)*size(p,3),1), 2)), ...
    real(sum( (reshape(permute(p,[1 3 2]),[],3)*jar_rotmat) .* repmat([0 -1 0],size(p,1)*size(p,3),1), 2)) ), size(p,1), 1, []) )<fid_mark_lengthphi) ) & ...
    ( (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(1))<fid_mark_pen) | ...
    (abs(sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2)-fid_mark_z(2))<fid_mark_pen) ) ) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'cylinder testmark';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1]*(jar_rotmat'), jar_cylrad-jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2) > 0) & ...
    (sum((p - repmat(test_cyl,[size(p,1),1,size(p,3)])).^2,2)<=testmark_radius^2), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'hemisphere testmark';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    jar_rotmat*[(jar_cylrad-jar_cylthick)^-2 0 0 ; 0 (jar_cylrad-jar_cylthick)^-2 0 ; 0 0 (jar_axrad-jar_axthick)^-2 ]*jar_rotmat', ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum(p.*repmat([0 0 1]*(jar_rotmat'),[size(p,1) 1 size(p,3)]),2) <= 0) & ...
    (sum((p - repmat(test_sph,[size(p,1),1,size(p,3)])).^2,2)<=testmark_radius^2), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%% create ray lists
[raydirections pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixel_pitch, .5*(1+cam_resolution), cam0_focallength, ...
    cam0_rotate(1), cam0_rotate(2), cam0_rotate(3), cam0_distortion*cam0_focallength^-2, lens_type);
rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{1} = pixelmap;
ray_startingpoints{1} = repmat(cam0_position,size(raydirections,1),1);

[raydirections pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixel_pitch, .5*(1+cam_resolution), cam1_focallength, ...
    cam1_rotate(1), cam1_rotate(2), cam1_rotate(3), cam1_distortion*cam1_focallength^-2, lens_type);
rays{2} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{2} = pixelmap;
ray_startingpoints{2} = repmat(cam1_position,size(raydirections,1),1);

%% find test points
% this will find the 4 ends and center of the top front fiducial mark, and
% the 2 horizontal ends and center of the bottom mark, in this order:
% [tc tt tb tl tr bc bl br]

testpoints = zeros(8,3);
testpoints(1,:) = jar_cylrad*[0, -1, 0] + [0 0 fid_mark_z(2)];
testpoints(2,:) = jar_cylrad*[0, -1, 0] + [0 0 fid_mark_z(2)+fid_mark_length];
testpoints(3,:) = jar_cylrad*[0, -1, 0] + [0 0 fid_mark_z(2)-fid_mark_length];
testpoints(4,:) = jar_cylrad*[-sin(fid_mark_lengthphi), -cos(fid_mark_lengthphi), 0] + [0 0 fid_mark_z(2)];
testpoints(5,:) = jar_cylrad*[sin(fid_mark_lengthphi), -cos(fid_mark_lengthphi), 0] + [0 0 fid_mark_z(2)];
testpoints(6,:) = jar_cylrad*[0, -1, 0] + [0 0 fid_mark_z(1)];
testpoints(7,:) = jar_cylrad*[-sin(fid_mark_lengthphi), -cos(fid_mark_lengthphi), 0] + [0 0 fid_mark_z(1)];
testpoints(8,:) = jar_cylrad*[sin(fid_mark_lengthphi), -cos(fid_mark_lengthphi), 0] + [0 0 fid_mark_z(1)];

testpoints = testpoints*(jar_rotmat');

for c=1:2
    test_pixels{c} = zeros(8,2);
    switch c
        case 1
            campos = cam0_position;
            camfoc = cam0_focallength;
            camrot = cam0_rotate;
            camdistort = cam0_distortion*cam0_focallength^-2;
        case 2
            campos = cam1_position;
            camfoc = cam1_focallength;
            camrot = cam1_rotate;
            camdistort = cam1_distortion*cam1_focallength^-2;
    end
    M1 = [cos(camrot(2)) -sin(camrot(2)) 0 ; sin(camrot(2)) cos(camrot(2)) 0 ; 0 0 1];
    M2 = [1 0 0 ; 0 cos(camrot(1)) -sin(camrot(1)) ; 0 sin(camrot(1)) cos(camrot(1))];
    M3 = [cos(camrot(3)) 0 sin(camrot(3)) ; 0 1 0 ; -sin(camrot(3)) 0 cos(camrot(3))];

    camrotmat = M1*M2*M3;
    
    for t=1:8
        xzplane = testpoints(t,[1 3]) - campos([1 3]);
        x = sqrt(sum((xzplane).^2));
        y = [window_inside-window_thickness-campos(2), window_thickness, testpoints(t,2)-window_inside];
        n = [n_air, n_glass, n_glycol];
        sintheta_guess = max(min(x / sum( ...
            y./(n/n(1))),.9),.1);
        opts = optimset('Jacobian','on','Display','off');
        sintheta = fsolve(@(sintheta)(AquariumSolver(sintheta,x,y,n)), sintheta_guess, opts);
        
        tantheta = sintheta/sqrt(1-sintheta^2);
        r = y(1)*tantheta;
        
        phi = atan2(xzplane(2),xzplane(1));
        
        p = [r*cos(phi), y(1), r*sin(phi)];
        p = p*camrotmat;
        p = -p * camfoc/p(2);
%         ray_direction = -pixel_location .* ...
%             (1 + radial_distortion .* repmat(sum(pixel_location(:,[1 3]).^2,2),1,3) .* ...
%             repmat([1 0 1],size(pixel_location,1),1));
        opts = optimset('Jacobian','off','Display','off');
        distortion_function = @(f)(f.*repmat(1+camdistort*sum(f.^2,2),1,2));
        p([1 3]) = fsolve(@(f)(distortion_function(f)-p([1 3])), p([1 3]), opts);
        test_pixels{c}(t,:) = (p([1 3])./cam_pixel_pitch) .* [-1 1] + .5*(1+cam_resolution);
        
    end
end

%% all done!
