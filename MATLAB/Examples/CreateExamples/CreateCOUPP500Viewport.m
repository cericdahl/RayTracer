
function [surface_list rays ray_startingpoints pixels] = CreateCOUPP500Viewport(geospecs)

%% set defaults
if nargin<1 || isempty(geospecs)
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

%% default dimensions
%%% indices of refraction
% n_hydraulic = 1.434;
% n_air = 1.00;
% n_pressurewindow = 1.52;
% 
% %%% viewport dims
% r0 = 12.5 * .5 * 2.54; % nozzle radius
% r1 = 12 * .5 * 2.54; % max radius of cone
% r2 = 6.5 * .5 * 2.54; % radius of cone at inside of window
% r3 = 2 * .5 * 2.54; % visible radius at outside of window
% a = pi/6; % half-angle of cone
% t = 2 * 2.54; % thickness of window
% 
% %%% viewport location
% s = 36 * 2.54; % cylinder radius to cone vertex
% z = 0; % height to cone vertex
% phi = 0; % azimuth to cone vertex
% theta = 0; % pitch of viewport (positive is looking upward)
% % viewport has no yaw
% 
% %%% camera location
% cam_setback = 1 * 2.54; % distance from window-outside to cam optical center
% % camera axis is viewport cone axis
% 
% %%% camera image details
% f = .8; % focal length of lens
% lens = 'theta'; % type of lens
% res = [1088 2048];
% ccd = [1088 2048] .* [5.5e-4 5.5e-4];
% camera pitch is theta;
% camera yaw is phi;

%% apply geospecs
fn = fieldnames(geospecs);
for n=1:length(fn)
    if ~isempty(geospecs.(fn{n}))
        eval([fn{n} '=geospecs.(fn{n});']);
    end
end

%% Make helpful bits
rotmat = [cos(phi), -sin(phi), 0 ; sin(phi), cos(phi), 0 ; 0, 0, 1] * ...
    [1, 0, 0, ; 0, cos(theta), -sin(theta) ; 0, sin(theta), cos(theta)];

cone_inward_axis = rotmat * [0 ; 1 ; 0];
cone_outward_axis = -cone_inward_axis;

cone_min_l = r2 / tan(a);
cone_max_l = r1 / tan(a);

cone_vertex = [s*sin(phi); -s*cos(phi); z];
cone_vertex_coneframe = rotmat' * cone_vertex;
cone_y = cone_vertex_coneframe(2);
cone_z = cone_vertex_coneframe(3);

cone_Q = rotmat * [1,0,0;0,-(tan(a)^2),0;0,0,1] * rotmat';
cone_P = rotmat * [0; 2*(tan(a)^2)*cone_y; -2*cone_z];
cone_R = cone_z^2 - (tan(a)*cone_y)^2;

vpname = sprintf('VP_%.0fdeg_%.0f"',180*phi/pi,z/2.54);

window_air_center = cone_vertex + cone_inward_axis*(cone_min_l-t);
window_glycol_center = cone_vertex + cone_inward_axis*cone_min_l;
nozzle_end_center = cone_vertex + cone_inward_axis*cone_max_l;

%% Make 7 surface

surface_list(end+1).description = [vpname, ' - nozzle'];
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    cone_vertex, cone_inward_axis, r0);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) >= (cone_min_l - t - cam_setback) ) & ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) <= cone_max_l ) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = [vpname, ' - tunnel'];
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    cone_vertex, cone_inward_axis, r2);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) >= (cone_min_l - t) ) & ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) <= cone_min_l ) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = [vpname, ' - baffle'];
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    cone_vertex, cone_inward_axis, r3);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) >= (cone_min_l - t - cam_setback) ) & ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) <= (cone_min_l - t) ) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = [vpname, ' - cone'];
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    cone_Q, cone_P, cone_R);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) >= cone_min_l ) & ...
    ( ((p(:,1,:) - cone_vertex(1))*cone_inward_axis(1) + (p(:,2,:) - cone_vertex(2))*cone_inward_axis(2) + (p(:,3,:) - cone_vertex(3))*cone_inward_axis(3)) <= cone_max_l ) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = [vpname, ' - glass-air interface'];
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    window_air_center, cone_outward_axis);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ( ( ((p(:,1,:)-window_air_center(1)).^2) + ((p(:,2,:)-window_air_center(2)).^2) + ((p(:,3,:)-window_air_center(3)).^2) ) <= (r2^2) ) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_pressurewindow;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = [vpname, ' - glycol-glass interface'];
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    window_glycol_center, cone_outward_axis);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ( ( ((p(:,1,:)-window_glycol_center(1)).^2) + ((p(:,2,:)-window_glycol_center(2)).^2) + ((p(:,3,:)-window_glycol_center(3)).^2) ) <= (r2^2) ) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = [vpname, ' - end-annulus'];
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    nozzle_end_center, cone_outward_axis);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ( ( ((p(:,1,:)-nozzle_end_center(1)).^2) + ((p(:,2,:)-nozzle_end_center(2)).^2) + ((p(:,3,:)-nozzle_end_center(3)).^2) ) >= (r1^2) ) & ...
    ( ( ((p(:,1,:)-nozzle_end_center(1)).^2) + ((p(:,2,:)-nozzle_end_center(2)).^2) + ((p(:,3,:)-nozzle_end_center(3)).^2) ) <= (r0^2) ) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%% create ray lists
[raydirections pixelmap] = GenerateRaysFromCamera(res, ccd./res, .5*(1+res), f, ...
    theta, phi, 0, barrel_d, lens);

rays = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels = pixelmap;
ray_startingpoints = repmat(cone_vertex(:)' + (cone_min_l-t-cam_setback)*cone_inward_axis(:)',size(raydirections,1),1);
