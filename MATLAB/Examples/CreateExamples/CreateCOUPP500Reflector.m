
function surface_list = CreateCOUPP500Reflector(geospecs)

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
% n_hydraulic = 1.434;
% 
% r1 = 10; % radius at bottom
% r2 = 10; % radius at top
% z1 = 0; % z at bottom
% z2 = 10; % if z2==z1 then this is a plane and we ignore r2
% phi = 0; % phi = 0 centers reflector on +y axis
% half_width = pi/2; % radians
% axis_offset = 0; % positive offset moves reflector farther from jar

%% apply geospecs
fn = fieldnames(geospecs);
for n=1:length(fn)
    if ~isempty(geospecs.(fn{n}))
        eval([fn{n} '=geospecs.(fn{n});']);
    end
end

%% Make helpful bits
reflector_backdir = [-sin(phi), cos(phi), 0];
reflector_center = axis_offset*reflector_backdir + [0, 0, z1];
cos_halfwidth = cos(half_width);

if z1~=z2
    rdiff_over_zdiff = (r2-r1)/(z2-z1);
else
    rdiff_over_zdiff = 0;
end

x0 = reflector_center(1);
y0 = reflector_center(2);

Q = [1, 0, 0 ; 0, 1, 0 ; 0, 0, -(rdiff_over_zdiff^2)];
P = [-2*x0 ; -2*y0 ; -2*rdiff_over_zdiff*( r1 - z1*rdiff_over_zdiff )];
R = x0^2 + y0^2 - (( r1 - z1*rdiff_over_zdiff )^2);

reflector_name = sprintf('Reflector %.0f-%.0f, %.0fdeg',z1,z2,phi);

%% Make surface

if z1~=z2
    surface_list(end+1).description = reflector_name;
    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
        Q, P, R);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        (p(:,3,:) >= z1) & ...
        (p(:,3,:) <= z2) & ...
        ( ( (p(:,1,:) - reflector_center(1))*reflector_backdir(1) + (p(:,2,:) - reflector_center(2))*reflector_backdir(2) ) >= ...
        (cos_halfwidth * (r1 + (p(:,3,:)-z1)*rdiff_over_zdiff)) ) ), ...
        size(p,1), []));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;
else
    surface_list(end+1).description = reflector_name;
    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
        reflector_center, [0 0 -1]);
    surface_list(end).inbounds_function = @(p)(reshape( ( ...
        ( ( (p(:,1,:) - reflector_center(1)).^2 + (p(:,2,:) - reflector_center(2)).^2 ) <= (r1^2)) & ...
        ( ( (p(:,1,:) - reflector_center(1))*reflector_backdir(1) + (p(:,2,:) - reflector_center(2))*reflector_backdir(2) ) >= (cos_halfwidth*r1) ) ), ...
        size(p,1), []));
    surface_list(end).n_outside = inf;
    surface_list(end).n_inside = n_hydraulic;
    surface_list(end).surface_type = 'retro';
    surface_list(end).absorption = 1;
end


