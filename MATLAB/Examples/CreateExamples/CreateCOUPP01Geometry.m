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

function [surface_list rays ray_startingpoints pixels] = CreateCOUPP01Geometry(geospecs)

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
n_CF3I = 1.20; % C3F8, actually
n_H2O = 1.33;
n_quartz = 1.458;
n_glycol = 1.33; % bath water, actually
n_air = 1.00;
n_glass = 1.491;

% jar dimensions in cm
jar_cylthick = .5*2.54*(1.69-.92); % thickness of cylinder wall
jar_axthick = .5*2.54*(1.69-.92); % thickness of sphere wall at apex
jar_cylrad = .5*2.54*1.69; % outer radius of cylinder
jar_axrad = .5*2.54*1.69; % outer radius of sphere (along cylinder axis)

% bath dimensions in cm
bath_bottom = -100;%-2.54*(.5*1.69+2.5);
bath_top = 100;%25;
bath_rad = 4.2*2.54;
plexi_thickness = .3*2.54;

% cf3i volume
cf3i_density = 1.38; % actually c3f8
cf3i_mass = 30;

% cameras, position relative to air-side center of viewport, +y towards
% chamber
cam_x = 0;
cam_y = -(5.55+4.2+.3)*2.54;
cam_z = (1.1-.5*1.69)*2.54;
cam_f = 1.2;
cam_barreld = 0;
cam_lenstype = 'theta';
cam_sensorsize = [491 656]*.00099;
cam_resolution = [491 656];

cam_pitch = 0;
cam_yaw = 0;
cam_roll = 0;

% % grid on reflector
% grid_xphase = 0;
% grid_zphase = 0;
% grid_minorlinehalfwidth = .05;
% grid_majorlinehalfwidth = .1;
% grid_majorpitch = 2.54;
% grid_minordivs = 4;

%% apply geospecs
fn = fieldnames(geospecs);
for n=1:length(fn)
    if ~isempty(geospecs.(fn{n}))
        eval([fn{n} '=geospecs.(fn{n});']);
    end
end

%% derived dimensions
liquid_level = ((cf3i_mass/cf3i_density) - (2/3)*pi*(jar_axrad-jar_axthick)*(jar_cylrad-jar_cylthick)^2) / ...
    (pi * (jar_cylrad-jar_cylthick)^2);

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

surface_list(end+1).description = 'outside surface of quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=0) & (p(:,3,:)<bath_top), size(p,1), [] ));
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

% surface_list(end+1).description = 'Bath bottom';
% surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
%     [0 0 bath_bottom], [0 0 -1]);
% surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) < (bath_rad.^2), size(p,1), [] ));
% surface_list(end).n_outside = inf;
% surface_list(end).n_inside = n_glycol;
% surface_list(end).surface_type = 'normal';
% surface_list(end).absorption = 1;
% 
% 
% surface_list(end+1).description = 'Bath top';
% surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
%     [0 0 bath_top], [0 0 1]);
% surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) < (bath_rad.^2), size(p,1), [] ));
% surface_list(end).n_outside = inf;
% surface_list(end).n_inside = n_glycol;
% surface_list(end).surface_type = 'normal';
% surface_list(end).absorption = 1;

surface_list(end+1).description = 'Bath ID';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], bath_rad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=bath_bottom) & (p(:,3,:)<bath_top), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;


surface_list(end+1).description = 'Bath OD';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], bath_rad+plexi_thickness);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=bath_bottom) & (p(:,3,:)<bath_top), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;


%% create ray lists
[raydirections pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixelpitch, .5*(1+cam_resolution), cam_f, ...
    cam_pitch*pi/180, cam_yaw*pi/180, cam_roll*pi/180, cam_barreld.*(cam_f.^(-2*(1:length(cam_barreld)))), cam_lenstype);

rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{1} = pixelmap;
ray_startingpoints{1} = repmat([cam_x cam_y cam_z],size(raydirections,1),1);

