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

function surface_list = CreateCOUPP500Jar(geospecs)

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

%% indices of refraction and dimensions used below
% n_target = 1.31;
% n_buffer = 1.33;
% n_jar = 1.458;
% n_hydraulic = 1.434;
% 
% % jar dimensions in cm
% jar_cylthick = .25; % thickness of cylinder wall
% jar_axthick = .25; % thickness of sphere wall at apex
% jar_cylrad = 7.5; % .5*OD of cylinder
% jar_axrad = 7.5; % .5*OD of sphere (along cylinder axis)
% 
% jar_cyllength = 7.62;
% jar_flangerad = 7.5;
% 
% liquid_level = 10;


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

%% create 7 surfaces

surface_list(end+1).description = 'inside surface of quartz cylinder below water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad - jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=0) & (p(:,3,:)<=liquid_level), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz cylinder above water';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad - jar_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=liquid_level) & (p(:,3,:)<=jar_cyllength), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_buffer;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], jar_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=0) & (p(:,3,:)<=jar_cyllength), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_hemi_inside_Q, ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<=0), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    quartz_hemi_outside_Q, ...
    [0 0 0], -1);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<=0), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'CF3I - water interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 liquid_level], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= ((jar_cylrad-jar_cylthick)^2), size(p,1), [] ));
surface_list(end).n_outside = n_buffer;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'jar top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 jar_cyllength], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= (jar_flangerad^2), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

