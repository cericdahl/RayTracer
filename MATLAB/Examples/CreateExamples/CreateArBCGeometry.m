% function [surface_list rays ray_startingpoints pixels] = Create30LGeometry()
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

function [surface_list, rays, ray_startingpoints, pixels] = CreateArBCGeometry(geospecs)

%% set defaults
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

%% useful equations
% n0 = index of refraction at known density
% r = rho/rho0 (rho = density now, rho0 = density corresponding to n0)
% returns index of refraction at density rho
clausius_mossotti = @(n0, r)(sqrt(((1 + 2*r).*n0.*n0 + 2 - 2*r)./((1 - r).*n0.*n0 + 2 + r)));

%% indices of refraction and dimensions used below
n_target = 1.17;% n=1.224 for Ar @ 940nm, 90K
n_jar = 1.4512; % SiO2 @ 940nm
n_hydraulic = 1.22;% 1.21 @ 940nm, 146K ;  1.237 @ 7eV, 146K, another said 1.515;
n_pressurewindow = 1.7569; % Al2O3 @ 940nm
n_pressurewall = inf;
n_air = 1.00;

% jar dimensions in cm
ojar_thick = .25; % thickness of cylinder wall
ojar_cylrad = 7.5; % outer radius of cylinder
ojar_axrad = 15; % outer radius of sphere (along cylinder axis)
ojar_knucklerad = 2.5;
ojar_cyllength = 40;
ojar_elevation = 20;

ijar_thick = .25; % thickness of cylinder wall
ijar_cylrad = 6.5; % outer radius of cylinder
ijar_axrad = 13; % outer radius of sphere (along cylinder axis)
ijar_knucklerad = 2.5;
ijar_cyllength = 20;
ijar_elevation = 0;

vp_s = 10; % radial position of air-side center of viewport
vp_elev = 60; % vertical position of air-side center of viewport
vp_win_rad = 1.73*.5*2.54; % radius of glass (black on circumference)
vp_air_rad = 1.25*.5*2.54; % radius of air-side can (black on circumference)
vp_can_rad = 2*2.54;
vp_can_wall = .125*2.54;
vp_flange_rad = 3.375*2.54;
vp_nip_rad = 1.75*.5*2.54; % radius of hydraulic-side nipple (black on circumference)
vp_win_thick = .25*2.54;
vp_nip_top = .5;
vp_theta = 6*pi/180;
vp_can_OAL = 6*2.54;
vp_flange_thick = [1,1,1,1,1]*2.54*.5;

rd_rad = 12;% reflector-diffuser radius
rd_top = 100;
rd_bot = 0;
rdcone_top = 120;
rdcone_toprad = 16;
rdtopcone_apex = 150;
rdtopcone_rad = 10.5;
rdtopcone_bot = -20;
rdbotcone_apex = -15.2;
rdbotcone_rad = 10.5;
rdbotcone_bot = -20;
        
pv_bot = -20;
pv_top = +100;
pv_rad = 30;
pv_thick = 1;
pv_axrad = 15;

% camera position, orientiation relative to viewport, -z is towards
% chamber, +y is towards jar axis
% (up is cos(vp_theta)\hat{z} + sin(vp_theta)\hat{y} )
cam_x = 0;
cam_y = 0;
cam_z = 5;
cam_f = .8;
cam_barreld = 0;
cam_lenstype = 'theta';
cam_sensorsize = [.1 .1];
cam_resolution = [480 640];

cam_pitch = 0;
cam_yaw = 0;
cam_roll = 0;

%% apply geospecs
fn = fieldnames(geospecs);
for n=1:length(fn)
    if ~isempty(geospecs.(fn{n}))
        eval([fn{n} '=geospecs.(fn{n});']);
    end
end

%% derived dimensions

t_o = [0 ojar_thick];
t_i = [0 ijar_thick];

r1 = [ojar_cylrad - t_o, ijar_cylrad - t_i];
r2 = [ojar_knucklerad - t_o, ijar_knucklerad - t_i];
r3 = [ojar_axrad - t_o, ijar_axrad - t_i];

s = r3.*(r1-r2)./(r3-r2); % axis to knuckle-dome transition

z = r2 .* sqrt(1 - (s./r3).^2); %  equator to knuckle-dome transition

d = r3 .* z .* ((1./r3)-(1./r2)); % equator to dome sphere center

cam_pixelpitch = cam_sensorsize ./ cam_resolution;

vp_axis = [0, -sin(vp_theta), cos(vp_theta)];
vp_center = [0, -vp_s, vp_elev];


head_out_Q = [pv_rad.^-2, 0, 0 ; 0, pv_rad.^-2, 0 ; 0, 0, pv_axrad.^-2];
head_in_Q = [(pv_rad-pv_thick).^-2, 0, 0 ; 0, (pv_rad-pv_thick).^-2, 0 ; 0, 0, (pv_axrad-pv_thick).^-2];
head_out_P = [0, 0, -2*pv_top .* pv_axrad.^-2];
head_in_P = [0, 0, -2*pv_top .* (pv_axrad-pv_thick).^-2];
head_out_R = (pv_top / pv_axrad).^2 - 1;
head_in_R = (pv_top / (pv_axrad - pv_thick)).^2 - 1;

rd_cone_b = (rdcone_toprad - rd_rad) ./ (rdcone_top - rd_top);
rd_cone_z0 = rd_top - (rd_rad / rd_cone_b);
rd_cone_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -rd_cone_b.^2];
rd_cone_P = [0, 0, 2*rd_cone_b^2*rd_cone_z0];
rd_cone_R = -(rd_cone_b * rd_cone_z0)^2;

rd_stcone_b = (rdcone_toprad - rdtopcone_rad) ./ (rdtopcone_bot - rdcone_top);
rd_stcone_z0 = rdtopcone_bot + (rdtopcone_rad / rd_stcone_b);
rd_stcone_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -rd_stcone_b.^2];
rd_stcone_P = [0, 0, 2*rd_stcone_b^2*rd_stcone_z0];
rd_stcone_R = -(rd_stcone_b * rd_stcone_z0)^2;

rd_topcone_b = rdtopcone_rad ./ (rdtopcone_apex - rdtopcone_bot);
rd_topcone_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -rd_topcone_b.^2];
rd_topcone_P = [0, 0, 2*rd_topcone_b^2*rdtopcone_apex];
rd_topcone_R = -(rd_topcone_b * rdtopcone_apex)^2;

rd_botcone_b = rdbotcone_rad ./ (rdbotcone_apex - rdbotcone_bot);
rd_botcone_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -rd_botcone_b.^2];
rd_botcone_P = [0, 0, 2*rd_botcone_b^2*rdbotcone_apex];
rd_botcone_R = -(rd_botcone_b * rdbotcone_apex)^2;

%% surface list

%%% first four silica cylinders
surface_list(end+1).description = 'inside surface of inner quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], r1(4));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<ijar_elevation) & (p(:,3,:)>=(ijar_elevation-ijar_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of inner quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], r1(3));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<ijar_elevation) & (p(:,3,:)>=(ijar_elevation-ijar_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of outer quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], r1(2));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<ojar_elevation) & (p(:,3,:)>=(ojar_elevation-ojar_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of outer quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], r1(1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<ojar_elevation) & (p(:,3,:)>=(ojar_elevation-ojar_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%%% then four silica domes
surface_list(end+1).description = 'inside surface of inner quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 ijar_elevation+d(4)], r3(4));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(4)+ijar_elevation)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of inner quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 ijar_elevation+d(3)], r3(3));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(3)+ijar_elevation)), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of outer quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 ojar_elevation+d(2)], r3(2));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(2)+ojar_elevation)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of outer quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 ojar_elevation+d(1)], r3(1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(1)+ojar_elevation)), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%%% now four silica knuckles
surface_list(end+1).description = 'inside surface of inner quartz knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 ijar_elevation], [0 0 1], r1(4)-r2(4), r2(4));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>ijar_elevation & p(:,3,:)<=(z(4)+ijar_elevation) & ...
    (p(:,1,:).^2+p(:,2,:).^2)>((r1(4)-r2(4))^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of inner quartz knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 ijar_elevation], [0 0 1], r1(3)-r2(3), r2(3));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>ijar_elevation & p(:,3,:)<=(z(3)+ijar_elevation) & ...
    (p(:,1,:).^2+p(:,2,:).^2)>((r1(3)-r2(3))^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of outer quartz knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 ojar_elevation], [0 0 1], r1(2)-r2(2), r2(2));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>ojar_elevation & p(:,3,:)<=(z(2)+ojar_elevation) & ...
    (p(:,1,:).^2+p(:,2,:).^2)>((r1(2)-r2(2))^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of outer quartz knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 ojar_elevation], [0 0 1], r1(1)-r2(1), r2(1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>ojar_elevation & p(:,3,:)<=(z(1)+ojar_elevation) & ...
    (p(:,1,:).^2+p(:,2,:).^2)>((r1(1)-r2(1))^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;


%%% now the viewport
surface_list(end+1).description = 'sight glass wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    vp_center, vp_axis, vp_air_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > 0) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_nip_top + vp_flange_thick(2))), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'camera can inner wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    vp_center, vp_axis, vp_can_rad - vp_can_wall);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (vp_nip_top + vp_flange_thick(2))) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_can_OAL + vp_flange_thick(2) + vp_nip_top)), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'camera can outer wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    vp_center, vp_axis, vp_can_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3))) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_can_OAL + vp_nip_top + vp_flange_thick(2) - vp_flange_thick(4))), size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'flange outer edge';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    vp_center, vp_axis, vp_flange_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (-vp_flange_thick(1)+vp_nip_top)) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)))) | ...
    ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4))) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_nip_top + vp_flange_thick(2) + vp_can_OAL + vp_flange_thick(5)))), size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'window wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    vp_center, vp_axis, vp_win_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= 0) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (-vp_win_thick)), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'window retainer outer wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    vp_center, vp_axis, vp_win_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= vp_nip_top) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > 0), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'pressure vessel nipple wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    vp_center, vp_axis, vp_nip_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= vp_nip_top) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (-vp_flange_thick(1)+vp_nip_top)), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'air side of viewport';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center, vp_axis);
surface_list(end).inbounds_function = @(p)(reshape( sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) <= (vp_air_rad^2), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'hydraulic side of viewport';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center - vp_axis*vp_win_thick, vp_axis);
surface_list(end).inbounds_function = @(p)(reshape( sum((p - repmat(vp_center - vp_axis*vp_win_thick,size(p,1),1,size(p,3))).^2,2) <= (vp_win_rad^2), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'viewport retainer';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center, vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) > (vp_air_rad^2)) & ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) <= (vp_win_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'nipple bottom';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center - vp_axis*(vp_flange_thick(1)-vp_nip_top), vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    (sum((p - repmat(vp_center - vp_axis*(vp_flange_thick(1)-vp_nip_top),size(p,1),1,size(p,3))).^2,2) > (vp_nip_rad^2)) & ...
    (sum((p - repmat(vp_center - vp_axis*(vp_flange_thick(1)-vp_nip_top),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'nipple top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center + vp_axis*vp_nip_top, vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    (sum((p - repmat(vp_center + vp_axis*vp_nip_top,size(p,1),1,size(p,3))).^2,2) > (vp_win_rad^2)) & ...
    (sum((p - repmat(vp_center + vp_axis*vp_nip_top,size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'can bot';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2)), vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2)),size(p,1),1,size(p,3))).^2,2) > (vp_air_rad^2)) & ...
    (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'can bot_top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)), vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)),size(p,1),1,size(p,3))).^2,2) > (vp_can_rad^2)) & ...
    (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'can top_bot';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4)), vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4)),size(p,1),1,size(p,3))).^2,2) > (vp_can_rad^2)) & ...
    (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'can top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL), vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'can very top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL + vp_flange_thick(5)), vp_axis);
surface_list(end).inbounds_function = @(p)( reshape(...
    sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL + vp_flange_thick(5)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2), size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%%% now other black surfaces to trap rays
surface_list(end+1).description = 'reflector/diffuser';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], rd_rad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>rd_bot) & ...
    (p(:,3,:)<=rd_top)) ,  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'reflector/diffuser cone';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    rd_cone_Q, rd_cone_P, rd_cone_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>rd_top) & ...
    (p(:,3,:)<rdcone_top)) ,  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'reflector/diffuser strip cone';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    rd_stcone_Q, rd_stcone_P, rd_stcone_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>rdcone_top) & ...
    (p(:,3,:)<rdtopcone_bot) & ...
    ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) - ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)),2).^2)) > (vp_nip_rad^2))) ,  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'reflector/diffuser topcone';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    rd_topcone_Q, rd_topcone_P, rd_topcone_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>rdtopcone_bot) & ...
    (p(:,3,:)<rdtopcone_apex)) ,  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'reflector/diffuser botcone';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    rd_botcone_Q, rd_botcone_P, rd_botcone_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>rdbotcone_bot) & ...
    (p(:,3,:)<rdbotcone_apex)) ,  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - cylinder outer wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], pv_rad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>pv_bot) & ...
    (p(:,3,:)<pv_top)) ,  ...
    size(p,1), []));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - cylinder inner wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], pv_rad - pv_thick);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>pv_bot) & ...
    (p(:,3,:)<pv_top)) ,  ...
    size(p,1), []));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - outer top';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    head_out_Q, head_out_P, head_out_R);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) > pv_top) & ...
    ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) - ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)),2).^2)) > (vp_flange_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - inner top';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    head_in_Q, head_in_P, head_in_R);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) > pv_top) & ...
    ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) - ...
    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)),2).^2)) > (vp_flange_rad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - bot';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, pv_bot], [0, 0, -1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= (pv_rad^2), size(p,1), [] ));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%% create ray lists
[raydirections, pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixelpitch, .5*(1+cam_resolution), cam_f, ...
    cam_pitch + vp_theta - (pi/2), cam_yaw, cam_roll, cam_barreld, cam_lenstype);

rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{1} = pixelmap;
ray_startingpoints{1} = repmat(vp_center + [cam_x, 0, 0] + cam_z*vp_axis + cam_y*cross(vp_axis, [1,0,0]),size(raydirections,1),1);

