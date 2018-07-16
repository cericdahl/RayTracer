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

function [surface_list rays ray_startingpoints pixels] = CreateNew2LGeometry(geospecs)

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
n_target = 1.31;
n_buffer = 1.33;
n_jar = 1.458;
n_hydraulic = 1.434;
n_pressurewindow = 1.52;
n_pressurewall = inf;
n_bath = 1.33;
n_bathwall = 1.33;
n_air = 1.00;

% jar dimensions in cm
jar_cylthick = .25; % thickness of cylinder wall
jar_axthick = .25; % thickness of sphere wall at apex
jar_cylrad = 7.5; % outer radius of cylinder
jar_axrad = 7.5; % outer radius of sphere (along cylinder axis)

jar_cyllength = 7.62;
jar_axrad_top = 7.5;
jar_axthick_top = .25;
jar_bellowsrad = 6.25;

% cf3i volume
target_mass = 4048;
target_density = 2;

% pressure vessel body (all dimensions inside-dimensions)
pv_cylbottom = -1*2.54;
pv_cyllength = 3.5*2.54;
pv_cylrad = 8.625*2.54; % inner radius of pressure vessel
pv_cylthick = 2.54*.375;
pv_axrad_top = 8.625*2.54;
pv_axrad_bot = 8.625*2.54;
pv_portrad_top = 3.03*2.54;
pv_portrad_bot = 3.03*2.54;
pv_top = 50.5; % height relative to jar center, from CM
pv_bot = -30; % pipe length to flange surface
pv_absorption = 1;

% pressure vessel view ports
vp_outerrad = .5*6.625*2.54;
vp_innerrad = .5*4*2.54;
vp_winrad = 4*.5*2.54;
vp_conelength = 3*2.54;%(2.52-2.5*.5)*2.54;
vp_innerlength = .3*2.54;
vp_winthick = 0.9*2.54;
vp_totallength = 12*2.54;
vp_height = .75*2.54;
vp_phi = 45*pi/180;
vp_lightring_innerrad = 2.54;
vp_lightring_outerrad = 2*2.54;

% reflector wall
tworeflectors = true;
ref_offaxis = 14.1;%5*2.54;
ref_cylrad = 30;%8.25*2.54;
ref_slope_top = 1;%8.25*2.54;
ref_slope_bot = 1;%8.25*2.54;
ref_azwidth = pi/2;
ref_cyllength = 3;%3.5*2.54;
ref_cylbottom = vp_height-.5*ref_cyllength;%-1*2.54;
ref_toplength = 6*2.54;
ref_botlength = 8*2.54;
ref_slope_bot2 = 2;%8.25*2.54;
ref_bot2length = 8*2.54;

% bath
bath_cylrad = 2.54*120; % outer radius of bath
bath_cylthick = 2.54*.375;
bath_cylbottom = -2.54;
bath_cyllength = 2.54*12;

% cameras, position relative to air-side center of viewport, +y towards
% chamber
cam_x = 0;
cam_y = -5;
cam_z = 0;
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
liquid_level = ( (target_mass / target_density) - ((2/3)*pi*(jar_cylrad-jar_cylthick)^2*(jar_axrad-jar_axthick)) ) / ...
    (pi*(jar_cylrad-jar_cylthick)^2);

if liquid_level > jar_cyllength
    liquid_level = jar_cyllength;
    disp('warning, tried to put too much CF3I in jar in CreateNew2LGeometry');
end

vp_rotmat = [cos(vp_phi) sin(vp_phi) 0 ; -sin(vp_phi) cos(vp_phi) 0; 0 0 1];

quartz_hemi_inside_Q = [(jar_cylrad-jar_cylthick)^-2 0 0 ; 0 (jar_cylrad-jar_cylthick)^-2 0 ; 0 0 (jar_axrad-jar_axthick)^-2 ];
quartz_hemi_outside_Q = [(jar_cylrad)^-2 0 0 ; 0 (jar_cylrad)^-2 0 ; 0 0 (jar_axrad)^-2 ];

quartz_upperhemi_inside_Q = [(jar_cylrad-jar_cylthick)^-2, 0, 0 ; 0, (jar_cylrad-jar_cylthick)^-2, 0 ; 0, 0, (jar_axrad_top-jar_axthick_top)^-2];
quartz_upperhemi_inside_P = [0, 0, -2*jar_cyllength*(jar_axrad_top-jar_axthick_top)^-2];
quartz_upperhemi_inside_R = (jar_cyllength/(jar_axrad_top-jar_axthick_top))^2 - 1;

quartz_upperhemi_outside_Q = [jar_cylrad^-2, 0, 0 ; 0, jar_cylrad^-2, 0 ; 0, 0, jar_axrad_top^-2];
quartz_upperhemi_outside_P = [0, 0, -2*jar_cyllength*jar_axrad_top^-2];
quartz_upperhemi_outside_R = (jar_cyllength/jar_axrad_top)^2 - 1;

top_dome_Q = [pv_cylrad^-2, 0, 0 ; 0, pv_cylrad^-2, 0 ; 0, 0, pv_axrad_top^-2];
top_dome_P = [0, 0, -2*(pv_cylbottom+pv_cyllength)*pv_axrad_top^-2];
top_dome_R = ((pv_cylbottom+pv_cyllength)/pv_axrad_top)^2 - 1;

bottom_dome_Q = [pv_cylrad^-2, 0, 0 ; 0, pv_cylrad^-2, 0 ; 0, 0, pv_axrad_bot^-2];
bottom_dome_P = [0, 0, -2*pv_cylbottom*pv_axrad_bot^-2];
bottom_dome_R = (pv_cylbottom/pv_axrad_bot)^2 - 1;

pv_botreflector = pv_cylbottom - pv_axrad_bot*sqrt(1-(pv_portrad_bot/pv_cylrad)^2);

cone_a2 = ((vp_outerrad-vp_innerrad)/vp_conelength)^2;
cone_y0 = vp_winthick+vp_innerlength-vp_totallength-(vp_conelength*vp_innerrad/(vp_outerrad-vp_innerrad));
cone_Q = [1, 0, 0 ; 0, -cone_a2, 0 ; 0, 0, 1];
cone_P = [0, 2*cone_a2*cone_y0, -2*vp_height];
cone_R = vp_height^2 - cone_a2*cone_y0^2;

reflector1_top_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -ref_slope_top^-2];
reflector1_top_P = [0, 2*ref_offaxis, 2*(ref_slope_top^-2)*(ref_cylbottom+ref_cyllength+ref_cylrad*ref_slope_top)];
reflector1_top_R = ref_offaxis^2 - ((ref_cylbottom+ref_cyllength+ref_cylrad*ref_slope_top)/ref_slope_top)^2;

reflector1_bot_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -ref_slope_bot^-2];
reflector1_bot_P = [0, 2*ref_offaxis, 2*(ref_slope_bot^-2)*(ref_cylbottom-ref_cylrad*ref_slope_bot)];
reflector1_bot_R = ref_offaxis^2 - ((ref_cylbottom-ref_cylrad*ref_slope_bot)/ref_slope_bot)^2;

reflector2_top_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -ref_slope_top^-2];
reflector2_top_P = [-2*sin(vp_phi)*ref_offaxis, 2*cos(vp_phi)*ref_offaxis, 2*(ref_slope_top^-2)*(ref_cylbottom+ref_cyllength+ref_cylrad*ref_slope_top)];
reflector2_top_R = ref_offaxis^2 - ((ref_cylbottom+ref_cyllength+ref_cylrad*ref_slope_top)/ref_slope_top)^2;

reflector2_bot_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -ref_slope_bot^-2];
reflector2_bot_P = [-2*sin(vp_phi)*ref_offaxis, 2*cos(vp_phi)*ref_offaxis, 2*(ref_slope_bot^-2)*(ref_cylbottom-ref_cylrad*ref_slope_bot)];
reflector2_bot_R = ref_offaxis^2 - ((ref_cylbottom-ref_cylrad*ref_slope_bot)/ref_slope_bot)^2;

reflector0_top_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -ref_slope_top^-2];
reflector0_top_P = [-2*sin(.5*vp_phi)*ref_offaxis, 2*cos(.5*vp_phi)*ref_offaxis, 2*(ref_slope_top^-2)*(ref_cylbottom+ref_cyllength+ref_cylrad*ref_slope_top)];
reflector0_top_R = ref_offaxis^2 - ((ref_cylbottom+ref_cyllength+ref_cylrad*ref_slope_top)/ref_slope_top)^2;

reflector0_bot_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -ref_slope_bot^-2];
reflector0_bot_P = [-2*sin(.5*vp_phi)*ref_offaxis, 2*cos(.5*vp_phi)*ref_offaxis, 2*(ref_slope_bot^-2)*(ref_cylbottom-ref_cylrad*ref_slope_bot)];
reflector0_bot_R = ref_offaxis^2 - ((ref_cylbottom-ref_cylrad*ref_slope_bot)/ref_slope_bot)^2;

reflector0_bot2_Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -ref_slope_bot2^-2];
reflector0_bot2_P = [-2*sin(.5*vp_phi)*ref_offaxis, 2*cos(.5*vp_phi)*ref_offaxis, 2*(ref_slope_bot2^-2)*(ref_cylbottom-ref_botlength-(ref_cylrad-(ref_botlength/ref_slope_bot))*ref_slope_bot2)];
reflector0_bot2_R = ref_offaxis^2 - ((ref_cylbottom-ref_botlength-(ref_cylrad-(ref_botlength/ref_slope_bot))*ref_slope_bot2)/ref_slope_bot2)^2;

cam_pixelpitch = cam_sensorsize ./ cam_resolution;

%% surface list

surface_list(end+1).description = 'inside surface of quartz cylinder below water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad - jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=0) & (p(:,3,:)<liquid_level), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz cylinder above water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad - jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=liquid_level) & (p(:,3,:)<jar_cyllength), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_buffer;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=0) & (p(:,3,:)<jar_cyllength), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_hemi_inside_Q, ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<0), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_hemi_outside_Q, ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<0), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'CF3I - water interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 liquid_level], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) < ((jar_cylrad-jar_cylthick)^2), size(p,1), [] ));
surface_list(end).n_outside = n_buffer;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of upper quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_upperhemi_inside_Q, quartz_upperhemi_inside_P, quartz_upperhemi_inside_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>=jar_cyllength) & ...
    ((p(:,1,:).^2+p(:,2,:).^2) >= (jar_bellowsrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_buffer;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of upper quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_upperhemi_outside_Q, quartz_upperhemi_outside_P, quartz_upperhemi_outside_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>=jar_cyllength) & ...
    ((p(:,1,:).^2+p(:,2,:).^2) >= (jar_bellowsrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Bellows cylinder (approx)';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_bellowsrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:) < pv_top) & ...
    (p(:,3,:) > jar_cyllength) & ...
    (( (p(:,1,:).^2 + p(:,2,:).^2)*(jar_cylrad-jar_cylthick)^-2 + (p(:,3,:)-jar_cyllength).^2*(jar_axrad_top-jar_axthick_top)^-2 ) >= 1 ) ), ...
    size(p,1), []));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - cylinder inside wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], pv_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>pv_cylbottom) & ...
    (p(:,3,:)<(pv_cylbottom+pv_cyllength)) & ...
    ((((p(:,3,:)-vp_height).^2 + p(:,1,:).^2) >= (vp_outerrad^2)) | (p(:,2,:)>0)) & ...
    ((((p(:,3,:)-vp_height).^2 + (p(:,1,:).*cos(vp_phi) + p(:,2,:).*sin(vp_phi)).^2) >= (vp_outerrad^2)) | ((p(:,2,:).*cos(vp_phi) - p(:,1,:).*sin(vp_phi))>0)) ),  ...
    size(p,1), []));
surface_list(end).n_outside = n_pressurewall;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = pv_absorption;

surface_list(end+1).description = 'PV - cylinder outside wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], pv_cylrad+pv_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>pv_bot) & ...
    (p(:,3,:)<pv_top) & ...
    ((((p(:,3,:)-vp_height).^2 + p(:,1,:).^2) >= (vp_outerrad^2)) | (p(:,2,:)>0)) & ...
    ((((p(:,3,:)-vp_height).^2 + (p(:,1,:).*cos(vp_phi) + p(:,2,:).*sin(vp_phi)).^2) >= (vp_outerrad^2)) | ((p(:,2,:).*cos(vp_phi) - p(:,1,:).*sin(vp_phi))>0)) ),  ...
    size(p,1), []));
surface_list(end).n_outside = n_bath;
surface_list(end).n_inside = n_pressurewall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = pv_absorption;

surface_list(end+1).description = 'PV - top dome';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    top_dome_Q, top_dome_P, top_dome_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:) >= (pv_cylbottom+pv_cyllength)) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2) >= (pv_portrad_top^2)) & ...
    ((((p(:,3,:)-vp_height).^2 + p(:,1,:).^2) >= (vp_outerrad^2)) | (p(:,2,:)>0)) & ...
    ((((p(:,3,:)-vp_height).^2 + (p(:,1,:).*cos(vp_phi) + p(:,2,:).*sin(vp_phi)).^2) >= (vp_outerrad^2)) | ((p(:,2,:).*cos(vp_phi) - p(:,1,:).*sin(vp_phi))>0)) ),  ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - top port pipe';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], pv_portrad_top);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:) < pv_top) & ...
    (p(:,3,:) > (pv_cylbottom+pv_cyllength)) & ...
    (( (p(:,1,:).^2 + p(:,2,:).^2)*pv_cylrad^-2 + (p(:,3,:)-(pv_cylbottom+pv_cyllength)).^2*pv_axrad_top^-2 ) >= 1 ) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'top flange';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 pv_top], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= (pv_portrad_top^2), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - bottom dome';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    bottom_dome_Q, bottom_dome_P, bottom_dome_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:) <= pv_cylbottom) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2) >= (pv_portrad_bot^2)) & ...
    ((((p(:,3,:)-vp_height).^2 + p(:,1,:).^2) >= (vp_outerrad^2)) | (p(:,2,:)>0)) & ...
    ((((p(:,3,:)-vp_height).^2 + (p(:,1,:).*cos(vp_phi) + p(:,2,:).*sin(vp_phi)).^2) >= (vp_outerrad^2)) | ((p(:,2,:).*cos(vp_phi) - p(:,1,:).*sin(vp_phi))>0)) ),  ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'PV - bottom port pipe';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], pv_portrad_bot);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:) > pv_bot) & ...
    (p(:,3,:) < pv_cylbottom) & ...
    (( (p(:,1,:).^2 + p(:,2,:).^2)*pv_cylrad^-2 + (p(:,3,:)-pv_cylbottom).^2*pv_axrad_bot^-2 ) >= 1 ) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'bottom flange';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 pv_bot], [0 0 -1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= (pv_portrad_bot^2), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'bath - cylinder inside wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], bath_cylrad-bath_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>bath_cylbottom) & ...
    (p(:,3,:)<(bath_cylbottom+bath_cyllength)) ),  ...
    size(p,1), []));
surface_list(end).n_outside = n_bathwall;
surface_list(end).n_inside = n_bath;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'bath - cylinder outside wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], bath_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>bath_cylbottom) & ...
    (p(:,3,:)<(bath_cylbottom+bath_cyllength)) ),  ...
    size(p,1), []));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_bathwall;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'VP1 - window casing';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 vp_height], [0 1 0], vp_winrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,2,:) > -vp_totallength) & ...
    (p(:,2,:) < (vp_winthick-vp_totallength)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP1 - lightring';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 -vp_totallength 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) > (vp_lightring_innerrad^2)) & ...
    ((p(:,1,:).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) <= (vp_lightring_outerrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP1 - glass-air interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 -vp_totallength 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2 + (p(:,3,:)-vp_height).^2) <= (vp_winrad^2)) & ...
    (((p(:,1,:).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) <= (vp_lightring_innerrad^2)) | ((p(:,1,:).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) > (vp_lightring_outerrad^2))) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'VP1 - glycol-glass interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 vp_winthick-vp_totallength 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2 + (p(:,3,:)-vp_height).^2) <= (vp_winrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'VP1 - end-annulus';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 vp_winthick-vp_totallength 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2 + (p(:,3,:)-vp_height).^2) > (vp_winrad^2)) & ...
    ((p(:,1,:).^2 + (p(:,3,:)-vp_height).^2) <= (vp_innerrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP1 - narrow pipe';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 vp_height], [0 1 0], vp_innerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,2,:) < (vp_winthick+vp_innerlength-vp_totallength)) & ...
    (p(:,2,:) > (vp_winthick-vp_totallength)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP1 - cone';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    cone_Q, cone_P, cone_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,2,:) >= (vp_winthick+vp_innerlength-vp_totallength)) & ...
    (p(:,2,:) <= (vp_winthick+vp_innerlength+vp_conelength-vp_totallength)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP1 - wide pipe';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 vp_height], [0 1 0], vp_outerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,2,:) > (vp_winthick+vp_innerlength+vp_conelength-vp_totallength)) & ...
    (p(:,2,:) < 0) & ...
    (((p(:,1,:).^2 + p(:,2,:).^2) > (pv_cylrad^2)) | (p(:,3,:) <= pv_cylbottom) | (p(:,3,:) >= (pv_cylbottom+pv_cyllength))) & ...
    ((( (p(:,1,:).^2 + p(:,2,:).^2)*pv_cylrad^-2 + (p(:,3,:)-(pv_cylbottom+pv_cyllength)).^2*pv_axrad_top^-2 ) >= 1 ) | (p(:,3,:) < (pv_cylbottom+pv_cyllength))) & ...
    ((( (p(:,1,:).^2 + p(:,2,:).^2)*pv_cylrad^-2 + (p(:,3,:)-pv_cylbottom).^2*pv_axrad_bot^-2 ) >= 1 ) | (p(:,3,:) > pv_cylbottom)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP2 - window casing';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 vp_height], [-sin(vp_phi) cos(vp_phi) 0], vp_winrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) > -vp_totallength) & ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) < (vp_winthick-vp_totallength)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP2 - lightring';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [vp_totallength*sin(vp_phi) -vp_totallength*cos(vp_phi) 0], [sin(vp_phi) -cos(vp_phi) 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) > (vp_lightring_innerrad^2)) & ...
    (((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) <= (vp_lightring_outerrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP2 - glass-air interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [vp_totallength*sin(vp_phi) -vp_totallength*cos(vp_phi) 0], [sin(vp_phi) -cos(vp_phi) 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-vp_height).^2) <= (vp_winrad^2)) & ...
    ((((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) <= (vp_lightring_innerrad^2)) | (((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-(vp_height+cam_z)).^2) > (vp_lightring_outerrad^2))) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'VP2 - glycol-glass interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [(vp_totallength-vp_winthick)*sin(vp_phi) (vp_winthick-vp_totallength)*cos(vp_phi) 0], [sin(vp_phi) -cos(vp_phi) 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-vp_height).^2) <= (vp_winrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'VP2 - end-annulus';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [(vp_totallength-vp_winthick)*sin(vp_phi) (vp_winthick-vp_totallength)*cos(vp_phi) 0], [sin(vp_phi) -cos(vp_phi) 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-vp_height).^2) > (vp_winrad^2)) & ...
    (((p(:,1,:)*cos(vp_phi) + p(:,2,:)*sin(vp_phi)).^2 + (p(:,3,:)-vp_height).^2) <= (vp_innerrad^2)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP2 - narrow pipe';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 vp_height], [-sin(vp_phi) cos(vp_phi) 0], vp_innerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) < (vp_winthick+vp_innerlength-vp_totallength)) & ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) > (vp_winthick-vp_totallength)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP2 - cone';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    vp_rotmat' * cone_Q * vp_rotmat, cone_P * vp_rotmat, cone_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) >= (vp_winthick+vp_innerlength-vp_totallength)) & ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) <= (vp_winthick+vp_innerlength+vp_conelength-vp_totallength)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'VP2 - wide pipe';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 vp_height], [-sin(vp_phi) cos(vp_phi) 0], vp_outerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) > (vp_winthick+vp_innerlength+vp_conelength-vp_totallength)) & ...
    ((p(:,2,:)*cos(vp_phi) - p(:,1,:)*sin(vp_phi)) < 0) & ...
    (((p(:,1,:).^2 + p(:,2,:).^2) > (pv_cylrad^2)) | (p(:,3,:) <= pv_cylbottom) | (p(:,3,:) >= (pv_cylbottom+pv_cyllength))) & ...
    ((( (p(:,1,:).^2 + p(:,2,:).^2)*pv_cylrad^-2 + (p(:,3,:)-(pv_cylbottom+pv_cyllength)).^2*pv_axrad_top^-2 ) >= 1 ) | (p(:,3,:) < (pv_cylbottom+pv_cyllength))) & ...
    ((( (p(:,1,:).^2 + p(:,2,:).^2)*pv_cylrad^-2 + (p(:,3,:)-pv_cylbottom).^2*pv_axrad_bot^-2 ) >= 1 ) | (p(:,3,:) > pv_cylbottom)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'bottom reflector';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 pv_botreflector], [0 0 -1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= (pv_portrad_bot^2), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;

if tworeflectors

    surface_list(end+1).description = 'reflector1 cylinder';
    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
        [0 -1 0]*ref_offaxis, [0 0 1], ref_cylrad);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)>=ref_cylbottom) & ...
        (p(:,3,:)<(ref_cylbottom+ref_cyllength)) & ...
        (((p(:,2,:)+ref_offaxis)./sqrt(p(:,1,:).^2+(p(:,2,:)+ref_offaxis).^2)) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector1 top cone';
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        reflector1_top_Q, reflector1_top_P, reflector1_top_R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)>=(ref_cylbottom+ref_cyllength)) & ...
        (p(:,3,:)<(ref_cylbottom+ref_cyllength+ref_toplength)) & ...
        (((p(:,2,:)+ref_offaxis)./sqrt(p(:,1,:).^2+(p(:,2,:)+ref_offaxis).^2)) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector1 bottom cone';
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        reflector1_bot_Q, reflector1_bot_P, reflector1_bot_R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)<ref_cylbottom) & ...
        (p(:,3,:)>(ref_cylbottom-ref_botlength)) & ...
        (((p(:,2,:)+ref_offaxis)./sqrt(p(:,1,:).^2+(p(:,2,:)+ref_offaxis).^2)) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector2 cylinder';
    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
        [sin(vp_phi) -cos(vp_phi) 0]*ref_offaxis, [0 0 1], ref_cylrad);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)>=ref_cylbottom) & ...
        (p(:,3,:)<(ref_cylbottom+ref_cyllength)) & ...
        ((cos(vp_phi)*((p(:,2,:)+cos(vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(vp_phi)*ref_offaxis).^2)) - sin(vp_phi)*((p(:,1,:)-sin(vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(vp_phi)*ref_offaxis).^2))) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector2 top cone';
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        reflector2_top_Q, reflector2_top_P, reflector2_top_R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)>=(ref_cylbottom+ref_cyllength)) & ...
        (p(:,3,:)<(ref_cylbottom+ref_cyllength+ref_toplength)) & ...
        ((cos(vp_phi)*((p(:,2,:)+cos(vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(vp_phi)*ref_offaxis).^2)) - sin(vp_phi)*((p(:,1,:)-sin(vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(vp_phi)*ref_offaxis).^2))) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector2 bottom cone';
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        reflector2_bot_Q, reflector2_bot_P, reflector2_bot_R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)<ref_cylbottom) & ...
        (p(:,3,:)>(ref_cylbottom-ref_botlength)) & ...
        ((cos(vp_phi)*((p(:,2,:)+cos(vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(vp_phi)*ref_offaxis).^2)) - sin(vp_phi)*((p(:,1,:)-sin(vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(vp_phi)*ref_offaxis).^2))) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

else
    
    surface_list(end+1).description = 'reflector0 cylinder';
    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
        [sin(.5*vp_phi) -cos(.5*vp_phi) 0]*ref_offaxis, [0 0 1], ref_cylrad);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)>=ref_cylbottom) & ...
        (p(:,3,:)<(ref_cylbottom+ref_cyllength)) & ...
        ((cos(.5*vp_phi)*((p(:,2,:)+cos(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2)) - sin(.5*vp_phi)*((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2))) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector0 top cone';
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        reflector0_top_Q, reflector0_top_P, reflector0_top_R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)>=(ref_cylbottom+ref_cyllength)) & ...
        (p(:,3,:)<(ref_cylbottom+ref_cyllength+ref_toplength)) & ...
        ((cos(.5*vp_phi)*((p(:,2,:)+cos(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2)) - sin(.5*vp_phi)*((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2))) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector0 first bottom cone';
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        reflector0_bot_Q, reflector0_bot_P, reflector0_bot_R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)<ref_cylbottom) & ...
        (p(:,3,:)>(ref_cylbottom-ref_botlength)) & ...
        ((cos(.5*vp_phi)*((p(:,2,:)+cos(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2)) - sin(.5*vp_phi)*((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2))) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;
    
    surface_list(end+1).description = 'reflector0 second bottom cone';
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        reflector0_bot2_Q, reflector0_bot2_P, reflector0_bot2_R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:)<ref_cylbottom-ref_botlength) & ...
        (p(:,3,:)>(ref_cylbottom-ref_botlength-ref_bot2length)) & ...
        ((cos(.5*vp_phi)*((p(:,2,:)+cos(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2)) - sin(.5*vp_phi)*((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis)./sqrt((p(:,1,:)-sin(.5*vp_phi)*ref_offaxis).^2+(p(:,2,:)+cos(.5*vp_phi)*ref_offaxis).^2))) > cos(.5*ref_azwidth)) ), ...
        size(p,1), [] ));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;
    
    surface_list(end+1).description = 'reflector_dummy';
    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
        [0 0 0], [0 0 1]);
    surface_list(end).inbounds_function = @(p)(false(size(p,1),size(p,3)));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

    surface_list(end+1).description = 'reflector_dummy';
    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
        [0 0 0], [0 0 1]);
    surface_list(end).inbounds_function = @(p)(false(size(p,1),size(p,3)));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;

end    
%% create ray lists
[raydirections pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixelpitch, .5*(1+cam_resolution), cam_f, ...
    cam_pitch, cam_yaw, cam_roll, cam_barreld, cam_lenstype);

rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{1} = pixelmap;
ray_startingpoints{1} = repmat([cam_x cam_y-vp_totallength cam_z+vp_height],size(raydirections,1),1);

