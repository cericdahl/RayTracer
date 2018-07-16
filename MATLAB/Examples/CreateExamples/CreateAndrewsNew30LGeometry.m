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

function [surface_list rays ray_startingpoints pixels] = CreateAndrewsNew30LGeometry()

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
n_glycol = 1.434;
n_air = 1.00;
n_glass = 1.52;
% n_lens = 1.52;

% all dimensions in cm
inner_radius = 14.6;
quartz_thickness = .4;
port_offset = (8)*2.54; % 0 is the original NuMI/D0 ports, positive offset = lower ports
cylinder_bottom = -19.75 + port_offset;
cylinder_top = cylinder_bottom-inner_radius-quartz_thickness+100.1-2.54-3.6;
cf3i_mass = 50000;
cf3i_volume = cf3i_mass / 2.096;
liquid_level = cylinder_bottom + ...
    ((cf3i_volume - (2*pi*(inner_radius^3)/3)) / (pi*(inner_radius^2)));
% vessel_radius = 12*2.54;
vessel_radius = (12-1.2)*2.54;
window_radius = 2*2.54;
tunnel_outer_radius = 3.99*2.54;
window_inside = -14*2.54 - .254;% - 1*2.54;
window_thickness = .9*2.54;
chamfer_vertex = -15.81*2.54;% - 1*2.54;
reflector_bottom = cylinder_bottom-14;
reflector_bottom_minor_radius = 9.4;
reflector_top = cylinder_top + 1;

tunnel_phi = (2*asin(tunnel_outer_radius / vessel_radius)) + (4/vessel_radius);

tunnel_rotmat = [cos(tunnel_phi) sin(tunnel_phi) 0 ; -sin(tunnel_phi) cos(tunnel_phi) 0; 0 0 1];

% lens_x = 1.55*2.54;
% lens_z = .42*2.54;
% lens_cylrad = .5*3.9;
% lens_minthick = .3;
% lens_rad = 2.8;

% fiber_lens_z = [3 -1];
% fiber_lens_cylrad = 1.5;
% fiber_lens_minthick = .3;
% fiber_lens_rad = 2.8;
% fiber_cylrad = 1;
% fiber_y = window_inside - window_thickness - fiber_lens_minthick - .5;

cam_x = 0;
cam_z = 0;
cam_y = window_inside - window_thickness - 3.5*2.54;
cam_focal_length = .65;
% cam_pixel_pitch = .0056265;
cam_ccd_dims = [.00055 .00055] .* ([1080 1080]-1);
cam_resolution = [108 108];
cam_pixel_pitch = cam_ccd_dims ./ (cam_resolution - 1);
% cam_distortion = [.23 .09 .23];

%% surface list

surface_list(end+1).description = 'inside surface of quartz cylinder below water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], inner_radius);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>cylinder_bottom) & (p(:,3,:)<liquid_level), size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz cylinder above water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], inner_radius);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=liquid_level) & (p(:,3,:)<cylinder_top), size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_H2O;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], inner_radius + quartz_thickness);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>cylinder_bottom) & (p(:,3,:)<cylinder_top), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 cylinder_bottom], inner_radius);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<=cylinder_bottom), size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 cylinder_bottom], inner_radius + quartz_thickness);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<=cylinder_bottom), size(p,1), [] ));
surface_list(end).n_outside = n_glycol;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'CF3I - water interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 liquid_level], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) < (inner_radius^2), size(p,1), [] ));
surface_list(end).n_outside = n_H2O;
surface_list(end).n_inside = n_CF3I;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'cylinder lid';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 cylinder_top], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= ((inner_radius+quartz_thickness)^2), size(p,1), [] ));
surface_list(end).n_outside = n_H2O;
surface_list(end).n_inside = n_H2O;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'chamfer';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    [1 0 0 ; 0 -1 0 ; 0 0 1], [0 2*chamfer_vertex 0], -chamfer_vertex^2);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2 + p(:,3,:).^2)>(window_radius^2)) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2)>(vessel_radius^2)) & ...
    ((p(:,1,:).^2 + p(:,3,:).^2)<(tunnel_outer_radius^2)) & ...
    (p(:,2,:) > chamfer_vertex), size(p,1), [] ));
surface_list(end).n_inside = n_glycol;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'tunnel';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 1 0], window_radius);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,2,:)>window_inside) & ...
    (p(:,2,:)<=(window_radius+chamfer_vertex)), size(p,1), [] ));
surface_list(end).n_inside = n_glycol;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'glass - glycol interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 window_inside 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2 + p(:,3,:).^2)<=(window_radius^2)), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'glass - air interface (plane)';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 window_inside-window_thickness 0], [0 -1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2 + p(:,3,:).^2)<=(window_radius^2)), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'reflector top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 reflector_top], [0 0 1]);
surface_list(end).inbounds_function = @(p)(true(size(p,1),size(p,3)));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'reflector bottom';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    [vessel_radius^-2, 0, 0 ; 0, vessel_radius^-2, 0 ; 0, 0, reflector_bottom_minor_radius^-2], ...
    [0 ; 0 ; -2*reflector_bottom*reflector_bottom_minor_radius^-2], (reflector_bottom/reflector_bottom_minor_radius)^2 - 1);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<=reflector_bottom), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'reflector wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], vessel_radius);
surface_list(end).inbounds_function = @(p)(reshape( ( (p(:,2,:)>0) | ...
    ((p(:,1,:).^2 + p(:,3,:).^2) > (tunnel_outer_radius^2)) ) & ...
    ( ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi))>0) | ...
    (((p(:,1,:)*cos(tunnel_phi) + p(:,2,:)*sin(tunnel_phi)).^2 + p(:,3,:).^2) > (tunnel_outer_radius^2)) ) & ...
    (p(:,3,:) > reflector_bottom) & (p(:,3,:) < reflector_top), ...
    size(p,1), [] ));
surface_list(end).n_inside = n_glycol;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'glass tunnel';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 1 0], window_radius);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,2,:)<=window_inside) & ...
    (p(:,2,:)>(window_inside-window_thickness)), size(p,1), [] ));
surface_list(end).n_inside = n_glass;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'tunnel wide portion';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 1 0], tunnel_outer_radius);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2 + p(:,2,:).^2)>(vessel_radius^2)) & ...
    (p(:,2,:) > (chamfer_vertex + tunnel_outer_radius)) & ...
    (p(:,2,:) < 0), ...
    size(p,1), [] ));
surface_list(end).n_inside = n_glycol;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'chamfer 2';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    tunnel_rotmat' * [1 0 0 ; 0 -1 0 ; 0 0 1] * tunnel_rotmat, [0 2*chamfer_vertex 0] * tunnel_rotmat, -chamfer_vertex^2);
surface_list(end).inbounds_function = @(p)(reshape( (((p(:,1,:)*cos(tunnel_phi) + p(:,2,:)*sin(tunnel_phi)).^2 + p(:,3,:).^2)>(window_radius^2)) & ...
    ((p(:,1,:).^2 + p(:,2,:).^2)>(vessel_radius^2)) & ...
    (((p(:,1,:)*cos(tunnel_phi) + p(:,2,:)*sin(tunnel_phi)).^2 + p(:,3,:).^2)<(tunnel_outer_radius^2)) & ...
    ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi)) > chamfer_vertex), size(p,1), [] ));
surface_list(end).n_inside = n_glycol;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'tunnel 2';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 1 0] * tunnel_rotmat, window_radius);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi))>window_inside) & ...
    ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi))<=(window_radius+chamfer_vertex)), size(p,1), [] ));
surface_list(end).n_inside = n_glycol;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'glass - glycol interface 2';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 window_inside 0] * tunnel_rotmat, [0 -1 0] * tunnel_rotmat);
surface_list(end).inbounds_function = @(p)(reshape( (((p(:,1,:)*cos(tunnel_phi) + p(:,2,:)*sin(tunnel_phi)).^2 + p(:,3,:).^2)<=(window_radius^2)), size(p,1), [] ));
surface_list(end).n_outside = n_glass;
surface_list(end).n_inside = n_glycol;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'glass - air interface (plane) 2';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 window_inside-window_thickness 0] * tunnel_rotmat, [0 -1 0] * tunnel_rotmat);
surface_list(end).inbounds_function = @(p)(reshape( (((p(:,1,:)*cos(tunnel_phi) + p(:,2,:)*sin(tunnel_phi)).^2 + p(:,3,:).^2)<=(window_radius^2)), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_glass;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'glass tunnel 2';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 1 0] * tunnel_rotmat, window_radius);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi))<=window_inside) & ...
    ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi))>(window_inside-window_thickness)), size(p,1), [] ));
surface_list(end).n_inside = n_glass;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'tunnel wide portion 2';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 1 0] * tunnel_rotmat, tunnel_outer_radius);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2 + p(:,2,:).^2)>(vessel_radius^2)) & ...
    ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi)) > (chamfer_vertex + tunnel_outer_radius)) & ...
    ((p(:,2,:)*cos(tunnel_phi) - p(:,1,:)*sin(tunnel_phi)) < 0), ...
    size(p,1), [] ));
surface_list(end).n_inside = n_glycol;
surface_list(end).n_outside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;


%% create ray lists
[raydirections pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixel_pitch, .5*(1+cam_resolution), cam_focal_length); %, ...
%     0, 0, 0, .5*cam_focal_length^-2);
rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{1} = pixelmap;
ray_startingpoints{1} = repmat([-cam_x cam_y cam_z],size(raydirections,1),1);

