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

function [surface_list rays ray_startingpoints pixels] = CreateCirteGeometry(geospecs)

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

%% indices of refraction and dimensions used below
n_CF3I = 1.31;
n_H2O = 1.33;
n_quartz = 1.458;
n_glycol = 1.33; % bath water, actually
n_air = 1.00;
n_glass = 1.491;

% jar dimensions in cm
jar_cylthick = .1; % thickness of cylinder wall
jar_axthick = .1; % thickness of sphere wall at apex
jar_cylrad = .5; % outer radius of cylinder
jar_axrad = .5; % outer radius of sphere (along cylinder axis)

% bath dimensions in cm
bath_bottom = -5;
bath_top = 25;
airgap_bottom = 5;
airgap_top = 7.5;
airgap_halfdepth = 5;
bath_halfwidth = 15;
bath_halfdepth = 15;

plexi_thickness = .3;

diffuser_standoff = 1;

% cf3i volume
liquid_level = 10;

% cameras, position relative to air-side center of viewport, +y towards
% chamber
cam_x = 0;
cam_y = -5;
cam_z = 0;
cam_f = .8;
cam_barreld = 0;
cam_lenstype = 'tan';
cam_sensorsize = [.4861 .6494];
cam_resolution = [491 656];

cam_pitch = 0;
cam_yaw = 0;
cam_roll = 0;

% grid on reflector
grid_xphase = 0;
grid_zphase = 0;
grid_minorlinehalfwidth = .05;
grid_majorlinehalfwidth = .1;
grid_majorpitch = 2.54;
grid_minordivs = 4;

%% apply geospecs
fn = fieldnames(geospecs);
for n=1:length(fn)
    if ~isempty(geospecs.(fn{n}))
        eval([fn{n} '=geospecs.(fn{n});']);
    end
end

%% derived dimensions

quartz_hemi_inside_Q = [(jar_cylrad-jar_cylthick)^-2 0 0 ; 0 (jar_cylrad-jar_cylthick)^-2 0 ; 0 0 (jar_axrad-jar_axthick)^-2 ];
quartz_hemi_outside_Q = [(jar_cylrad)^-2 0 0 ; 0 (jar_cylrad)^-2 0 ; 0 0 (jar_axrad)^-2 ];

cam_pixelpitch = cam_sensorsize ./ cam_resolution;

%% surface list

surface_list(end+1).description = 'inside surface of quartz cylinder below water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad - jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=0) & (p(:,3,:)<liquid_level), size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz cylinder above water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad - jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=liquid_level) & (p(:,3,:)<bath_top), size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_H2O;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder, bathlow';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=0) & (p(:,3,:)<(airgap_bottom-plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder, plexilow';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=(airgap_bottom-plexi_thickness)) & (p(:,3,:)<airgap_bottom), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder, airgap';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=airgap_bottom) & (p(:,3,:)<airgap_top), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder, plexihigh';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=airgap_top) & (p(:,3,:)<(airgap_top+plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder, bathhigh';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=(airgap_top+plexi_thickness)) & (p(:,3,:)<bath_top), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_hemi_inside_Q, ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<0), size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_hemi_outside_Q, ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<0), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'CF3I - water interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 liquid_level], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) < ((jar_cylrad-jar_cylthick)^2), size(p,1), [] ));
surface_list(end).n_outside = n_H2O;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Bath bottom';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 bath_bottom], [0 0 -1]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < bath_halfwidth ) & (abs(p(:,2,:)) < bath_halfdepth), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Plexi bottom';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 bath_bottom-plexi_thickness], [0 0 -1]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness) ) & (abs(p(:,2,:)) < (bath_halfdepth + plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'Bath top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 bath_top], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness) ) & (abs(p(:,2,:)) < (bath_halfdepth + plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'Bath front';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 -bath_halfdepth 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < bath_halfwidth ) & (p(:,3,:) < bath_top) & (p(:,3,:) > bath_bottom), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Plexi front';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 -bath_halfdepth-plexi_thickness 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & (p(:,3,:) < bath_top) & (p(:,3,:) > (bath_bottom - plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Bath back';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 bath_halfdepth 0], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < bath_halfwidth ) & (p(:,3,:) < bath_top) & (p(:,3,:) > bath_bottom), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Plexi back';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 bath_halfdepth+plexi_thickness 0], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & (p(:,3,:) < bath_top) & (p(:,3,:) > (bath_bottom - plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Bath left';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [-bath_halfwidth 0 0], [-1 0 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,2,:)) < bath_halfdepth ) & (p(:,3,:) < bath_top) & (p(:,3,:) > bath_bottom) & ...
    ~( (abs(p(:,2,:)) < (airgap_halfdepth + plexi_thickness)) & (p(:,3,:) > (airgap_bottom - plexi_thickness)) & (p(:,3,:) < (airgap_top + plexi_thickness)) ), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Plexi left';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [-bath_halfwidth-plexi_thickness 0 0], [-1 0 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,2,:)) < (bath_halfdepth + plexi_thickness)) & (p(:,3,:) < bath_top) & (p(:,3,:) > (bath_bottom - plexi_thickness)) & ...
    ~( (abs(p(:,2,:)) < airgap_halfdepth) & (p(:,3,:) > airgap_bottom) & (p(:,3,:) < airgap_top) ), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Bath right';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [bath_halfwidth 0 0], [1 0 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,2,:)) < bath_halfdepth ) & (p(:,3,:) < bath_top) & (p(:,3,:) > bath_bottom) & ...
    ~( (abs(p(:,2,:)) < (airgap_halfdepth + plexi_thickness)) & (p(:,3,:) > (airgap_bottom - plexi_thickness)) & (p(:,3,:) < (airgap_top + plexi_thickness)) ), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Plexi right';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [bath_halfwidth+plexi_thickness 0 0], [1 0 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,2,:)) < (bath_halfdepth + plexi_thickness)) & (p(:,3,:) < bath_top) & (p(:,3,:) > (bath_bottom - plexi_thickness)) & ...
    ~( (abs(p(:,2,:)) < airgap_halfdepth) & (p(:,3,:) > airgap_bottom) & (p(:,3,:) < airgap_top) ), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-inside front';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 -airgap_halfdepth 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & (p(:,3,:) < airgap_top) & (p(:,3,:) > airgap_bottom), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-outside front';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 -airgap_halfdepth-plexi_thickness 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < bath_halfwidth) & (p(:,3,:) < (airgap_top+plexi_thickness)) & (p(:,3,:) > (airgap_bottom-plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-inside back';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 airgap_halfdepth 0], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & (p(:,3,:) < airgap_top) & (p(:,3,:) > airgap_bottom), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-outside back';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 airgap_halfdepth+plexi_thickness 0], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < bath_halfwidth) & (p(:,3,:) < (airgap_top+plexi_thickness)) & (p(:,3,:) > (airgap_bottom-plexi_thickness)), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-inside top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 airgap_top], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & (abs(p(:,2,:)) < airgap_halfdepth) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2) >= (jar_cylrad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-outside top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 airgap_top+plexi_thickness], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < bath_halfwidth) & (abs(p(:,2,:)) < (airgap_halfdepth + plexi_thickness)) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2) >= (jar_cylrad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-inside bottom';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 airgap_bottom], [0 0 -1]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & (abs(p(:,2,:)) < airgap_halfdepth) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2) >= (jar_cylrad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Tunnel-outside bottom';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 airgap_bottom-plexi_thickness], [0 0 -1]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < bath_halfwidth) & (abs(p(:,2,:)) < (airgap_halfdepth + plexi_thickness)) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2) >= (jar_cylrad^2)), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Back Plane';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 bath_halfdepth+plexi_thickness+diffuser_standoff 0], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & ...
    (p(:,3,:) < bath_top) & ...
    (p(:,3,:) > (bath_bottom - plexi_thickness)) & ...
    ~( (abs(mod(p(:,1,:) + grid_xphase + .5*grid_majorpitch, grid_majorpitch) - .5*grid_majorpitch) < grid_majorlinehalfwidth) | ...
    (abs(mod(p(:,3,:) + grid_zphase + .5*grid_majorpitch, grid_majorpitch) - .5*grid_majorpitch) < grid_majorlinehalfwidth) | ...
    (abs(mod(p(:,1,:) + grid_xphase + .5*(grid_majorpitch/grid_minordivs), grid_majorpitch/grid_minordivs) - .5*(grid_majorpitch/grid_minordivs)) < grid_minorlinehalfwidth) | ...
    (abs(mod(p(:,3,:) + grid_zphase + .5*(grid_majorpitch/grid_minordivs), grid_majorpitch/grid_minordivs) - .5*(grid_majorpitch/grid_minordivs)) < grid_minorlinehalfwidth) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'Back Plane w/ grid';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 bath_halfdepth+plexi_thickness+diffuser_standoff 0], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( (abs(p(:,1,:)) < (bath_halfwidth + plexi_thickness)) & ...
    (p(:,3,:) < bath_top) & ...
    (p(:,3,:) > (bath_bottom - plexi_thickness)) & ...
    ~( (abs(mod(p(:,1,:) + grid_xphase + .5*grid_majorpitch, grid_majorpitch) - .5*grid_majorpitch) < grid_majorlinehalfwidth) | ...
    (abs(mod(p(:,3,:) + grid_zphase + .5*grid_majorpitch, grid_majorpitch) - .5*grid_majorpitch) < grid_majorlinehalfwidth) | ...
    (abs(mod(p(:,1,:) + grid_xphase + .5*(grid_majorpitch/grid_minordivs), grid_majorpitch/grid_minordivs) - .5*(grid_majorpitch/grid_minordivs)) < grid_minorlinehalfwidth) | ...
    (abs(mod(p(:,3,:) + grid_zphase + .5*(grid_majorpitch/grid_minordivs), grid_majorpitch/grid_minordivs) - .5*(grid_majorpitch/grid_minordivs)) < grid_minorlinehalfwidth) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%% create ray lists
[raydirections pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixelpitch, .5*(1+cam_resolution), cam_f, ...
    cam_pitch, cam_yaw, cam_roll, cam_barreld.*(cam_f.^(-2*(1:length(cam_barreld)))), cam_lenstype);

rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{1} = pixelmap;
ray_startingpoints{1} = repmat([cam_x cam_y-bath_halfdepth-plexi_thickness cam_z],size(raydirections,1),1);

