% function ray_interfaces = RayTracer(ray_startingpoints, rays, surfacelist, max_scatters, min_travel_length, follow_threshold, tir_handling)
%
% RayTracer propagates a set of rays through a geometry of reflecting and
% refracting surfaces.  The function inputs a set if initial rays and a
% list of surfaces.  It loops over surfaces, finding the intersection point
% of each ray with each surface, to determine the next scattering point of
% each ray.  Reflected and refracted rays are generated according to the
% surface properties.  The function iterates with the new set of rays,
% until all rays drop below threshold or the maximum number of scatters is
% reached.  The history of each ray is output in a structure array.
%
% inputs:
%           ray_startingpoints  -  N-by-3 matrix, where N is the number of
%                                    initial rays to follow, giving the
%                                    starting point for each ray.
%           rays                -  N-by-10 matrix giving the initial
%                                    direction, intensity, and polarization
%                                    of each ray.  The first 3 columns give
%                                    the forward direction of the ray
%                                    (these will be normalized if they
%                                    aren't already), columns 4-6 give a
%                                    direction non-parallel to the ray that
%                                    defines the s1 polarization axis
%                                    (these will be made normal to the ray
%                                    direction and normalized, if they
%                                    aren't already), and columns 7-10 are
%                                    the stokes parameters s0-s3, giving
%                                    the intensity and polarization of the
%                                    ray ( s0 gives the total intensity,
%                                    and s0^2 >= s1^2 + s2^2 + s3^2, see
%                                    7.2 in Jackson for more details)
%           surfacelist         -  A structure array defining the geometry
%                                    of scattering surfaces -- see
%                                    Create2LGeometry for details
%           max_scatters        -  The maximum number of scatters to
%                                    propagate rays through (the simulation
%                                    may stop before this point if there
%                                    are no rays above threshold left to
%                                    follow).
%           min_travel_length   -  A minimum travel length between
%                                    scatters.  This prevents rounding
%                                    errors from causing a ray to scatter
%                                    multiple times at the same point.
%                                    Rays with legitimate travel lengths
%                                    below this value will be INCORRECTLY
%                                    RECONSTRUCTED, so keep this small
%                                    (~1e-5 times your smallest dimension
%                                    is probably sufficient for most
%                                    geometries.)
%           follow_threshold    -  Refracted or reflected rays with an s0
%                                    below follow_threshold(1) or
%                                    follow_threshold(2), respectively,
%                                    will not be followed.  If
%                                    follow_threshold is a scalar, the same
%                                    threshold is used for both.
%           tir_handling        -   This determines what the refracted_rays
%                                     output is in the case of total internal
%                                     reflection.  The default (-1) gives a
%                                     refracted ray tangent to the surface with
%                                     zero intensity.  Any value >=0 will give
%                                     a ray with the same direction and
%                                     polarization as the reflected ray, with
%                                     intensity equal to the reflected
%                                     intensity times tir_handling.  This lets
%                                     you treat tir-rays like refracted rays,
%                                     which can be handy in geometry sims.
%                                     NOTE -- if follow_threshold(2) is
%                                     bigger than max(rays(:,7)) then
%                                     default tir_handling=1.
%
% output:
%       ray_interfaces  -  a structure array, where the array index is the
%                            scatter number.  Each element in the array has
%                            the following fields:
%           incoming_ray        -  M-by-10 array, where M is the number of
%                                    rays scattering in this iteration,
%                                    giving the direction, intensity, and
%                                    polarization of the incoming rays.
%                                    Rays that do not scatter are not
%                                    reported (to report all rays, enclose
%                                    your geometry in an absorbing box, for
%                                    example).
%           reflected_ray       -  M-by-10 array, giving the direction,
%                                    intensity, and polarization of the
%                                    reflected rays
%           refracted_ray       -  M-by-10 array, giving the direction,
%                                    intensity, and polarization of the
%                                    refracted rays
%           intersection_point  -  M-by-3 array, giving the points where
%                                    the incoming rays scatter
%           surface_normal      -  M-by-3 array, giving the
%                                    backward-pointing surface normal at
%                                    the intersection_point
%           ray_index           -  M-by-1 vector, giving the index of the
%                                    incoming ray (the input rays are
%                                    numbered 1:N) -- a negative ray_index
%                                    means the ray has undergone at least
%                                    one reflection in its history, so
%                                    there will be at most one ray with a
%                                    given positive index
%           surface_index       -  M-by-1 vector, giving the index of the
%                                    scattering surface hit by each
%                                    incoming_ray, where the index
%                                    indicates an element of the input
%                                    surfacelist.  Negative values indicate
%                                    an outward ray, positive values an
%                                    inward ray, as defined in the
%                                    surface geometry.
%           distace_traveled    -  M-by-1 vector, giving the distance
%                                    traveled by the ray since its last
%                                    scatter
%           n_incident          -  M-by-1 vector, giving the index of
%                                    refraction for the incoming ray
%           n_transmitted       -  M-by-1 vector, giving the index of
%                                    refraction for the refracted ray
% 
% 
% 12/17/09, CED

function ray_interfaces = RayTracer(ray_startingpoints, rays, surfacelist, max_scatters, min_travel_length, follow_threshold, tir_handling)

ray_interfaces = struct( ...
    'incoming_ray', {}, ...
    'reflected_ray', {}, ...
    'refracted_ray', {}, ...
    'intersection_point', {}, ...
    'surface_normal', {}, ...
    'ray_index', {}, ...
    'surface_index', {}, ...
    'distance_traveled', {}, ...
    'n_incident', {}, ...
    'n_transmitted', {});

%% set defaults
if nargin<7 || isempty(tir_handling)
    tir_handling = [];
end
if nargin<6 || isempty(follow_threshold)
    follow_threshold = 0;
elseif numel(follow_threshold)==1
    follow_threshold = [0;0] + follow_threshold;
else
    follow_threshold = follow_threshold(1:2);
end

if nargin<5 || isempty(min_travel_length)
    min_travel_length = eps;
else
    min_travel_length = min_travel_length(1);
end

if nargin<4 || isempty(max_scatters)
    max_scatters = 10;
else
    max_scatters = max_scatters(1);
end

%% check inputs
if nargin<3 || ~isstruct(surfacelist) || length(size(ray_startingpoints))~=2 || ...
        length(size(rays))~=2 || size(ray_startingpoints,2)~=3 || ...
        size(rays,2)~=10 || size(ray_startingpoints,1)~=size(rays,1)
    disp('Impropper input to RayTracer');
    return
end

numrays = size(rays,1);
rays(:,1:3) = rays(:,1:3) ./ repmat(abs(sqrt(sum(rays(:,1:3).^2, 2))), 1, 3);
rays(:,4:6) = rays(:,4:6) ./ repmat(abs(sqrt(sum(rays(:,4:6).^2, 2))), 1, 3);

%% now really set default tir_handling
if isempty(tir_handling)
    if follow_threshold(2) > max(rays(:,7))
        tir_handling = 1;
    else
        tir_handling = -1;
    end
end

%% initialize raylist
p_start = ray_startingpoints;
incoming_rays = rays;
ray_index = (1:numrays)';

%%  follow rays
num_scatters = 0;
while ~isempty(ray_index)
    if num_scatters >= max_scatters
        break
    end
    num_scatters = num_scatters + 1;
    
    %% find next scattering surface for each ray
    % initialize the variables containing scatter properties
    p_next = zeros(size(p_start));            % scatter point
    l_next = zeros(size(p_start,1),1) + inf;  % distance to scatter
    s_next = zeros(size(p_start));            % optical surface normal at scatter
    sm_next = zeros(size(p_start));           % mechanical surface normal at scatter
    n_next = zeros(size(p_start,1),2);        % indices of refraction at scatter
    abs_next = zeros(size(p_start,1),1);      % absorption at scatter
    six_next = zeros(size(p_start,1),1);      % scatter surface index
    surfacetype_next = zeros(size(p_start,1),1); % surface type

    % loop over surfaces 
    for n=1:length(surfacelist)
        % find scatter point, normal, distance, and orientation for
        % scatters on this surface
        [p_intersect, s_normal, l_ray, s_orientation] = ...
            surfacelist(n).intersect_function(p_start,incoming_rays(:,1:3));
        
        % save the actual normal vector, because that will be important
        s_normal_mech = s_normal; %N x 3 x M matrix
        
        % adjust the surface normal as appropriate for the surface_type
        int_surfacetype = 0;
        switch surfacelist(n).surface_type
            case 'diffuse' % reflect incoming ray into 2pi
                int_surfacetype = 1;
                                                
            case 'retro' % reflect incoming ray back where it came from
                s_normal = repmat(-incoming_rays(:,1:3), [1 1 size(p_intersect,3)]);
        end
        
        % find the closest intersection point in the inbounds portion of
        % this surface (counting only positive, real distances and ignoring
        % glancing blows)
        valid_intersection = ...
            ( surfacelist(n).inbounds_function(p_intersect) ) & ...
            ( imag(l_ray)==0 ) & ...
            ( s_orientation ~= 0 ) & ...
            ( ~isnan(l_ray) ) & ...
            ( l_ray < inf ) & ...
            ( l_ray > min_travel_length );
        l_ray(~valid_intersection) = inf;
        [l_ray  ix] = min(l_ray, [], 2);
        
        % find the intersection points, surface normals, and orientations
        % assosciated with these scatters -- this is all just array
        % manipulation to extract the rows assosciated with the scatter
        % points found above
        ixlist = repmat(reshape(1:size(p_intersect,3),1,[],1),[3 1 size(p_intersect,1)]);
        ixcut = ixlist == repmat(reshape(ix,1,1,[]),[3 size(p_intersect,3) 1]);
        p_intersect = permute(p_intersect,[2 3 1]);
        p_intersect = reshape(p_intersect(ixcut),3,[])';
        s_normal = permute(s_normal,[2 3 1]);
        s_normal = reshape(s_normal(ixcut),3,[])';
        s_normal_mech = permute(s_normal_mech,[2 3 1]);
        s_normal_mech = reshape(s_normal_mech(ixcut),3,[])';

        ixlist = repmat(reshape(1:size(s_orientation,2),[],1),[1 size(s_orientation,1)]);
        ixcut = ixlist == repmat(ix(:)',[size(s_orientation,2) 1]);
        s_orientation = s_orientation';
        s_orientation = reshape(s_orientation(ixcut),[],1);
        
        n_before_after = repmat([surfacelist(n).n_outside surfacelist(n).n_inside],length(s_orientation),1);
        n_before_after(s_orientation<0,:) = repmat([surfacelist(n).n_inside surfacelist(n).n_outside],sum(s_orientation<0),1);
        
        % if this is the closest scatter so far, update the scatter
        % property variables
        scatter_here = l_ray < l_next;
        l_next(scatter_here) = l_ray(scatter_here);
        s_next(scatter_here,:) = s_normal(scatter_here,:);
        sm_next(scatter_here,:) = s_normal_mech(scatter_here,:);
        p_next(scatter_here,:) = p_intersect(scatter_here,:);
        n_next(scatter_here,:) = n_before_after(scatter_here,:);
        abs_next(scatter_here) = surfacelist(n).absorption;
        six_next(scatter_here) = n .* s_orientation(scatter_here);
        surfacetype_next(scatter_here) = int_surfacetype;
    end
    
    %% refigure s_normal for rays hitting weird surfaces like diffuse
    diffuse_cut = surfacetype_next==1;
    n_diffuse = sum(diffuse_cut);
    
    cos_theta = sqrt(rand(n_diffuse,1)); % diffuse ~= isotropic
    sin_theta = sqrt(1-cos_theta.^2);
    phi = rand(n_diffuse,1) * 2 * pi;
    
    x_tmp = cross(s_next(diffuse_cut,:), repmat([1 0 0],n_diffuse,1), 2);
    y_tmp = cross(s_next(diffuse_cut,:), repmat([0 1 0],n_diffuse,1), 2);
    tmpcut = all(x_tmp==0,2);
    x_tmp(tmpcut,:) = y_tmp(tmpcut,:);
    x_tmp = x_tmp ./ repmat(abs(sqrt(sum(x_tmp.^2, 2))), 1, 3);
    y_tmp = cross(s_next(diffuse_cut,:), x_tmp);
    
    outdir = s_next(diffuse_cut,:).*repmat(cos_theta,1,3) + ...
        x_tmp.*repmat(sin_theta .* cos(phi),1,3) + ...
        y_tmp.*repmat(sin_theta .* sin(phi),1,3);
    
    s_normal_tmp = outdir - incoming_rays(diffuse_cut,1:3);
    s_next(diffuse_cut,:) = s_normal_tmp ./ repmat(abs(sqrt(sum(s_normal_tmp.^2, 2))), 1, 3);

    %% reflect and refract at the surface
    % not all rays scatter -- we ignore (and do not report) those that don't
    scatter_cut = l_next < inf;
    if ~any(scatter_cut)
        ray_index = [];
        continue
    end
    [refracted_rays reflected_rays] = RefractionReflectionAtInterface(incoming_rays(scatter_cut,:), ...
        s_next(scatter_cut,:), n_next(scatter_cut,1), n_next(scatter_cut,2), tir_handling);
    % apply the absorption coefficient
    refracted_rays(:,7:10) = refracted_rays(:,7:10) .* repmat(1-abs_next(scatter_cut),1,4);
    reflected_rays(:,7:10) = reflected_rays(:,7:10) .* repmat(1-abs_next(scatter_cut),1,4);

    %% we need to reindex the variables used in constructing the new rays
    ray_index = ray_index(scatter_cut);
    p_next = p_next(scatter_cut,:);
    
    %% store output
    
    % we round the indices because the sign command used to find
    % orientation sometimes gives non-integer surface indices
    ray_interfaces(num_scatters).incoming_ray = incoming_rays(scatter_cut,:);
    ray_interfaces(num_scatters).refracted_ray = refracted_rays;
    ray_interfaces(num_scatters).reflected_ray = reflected_rays;
    ray_interfaces(num_scatters).intersection_point = p_next;
    ray_interfaces(num_scatters).surface_normal = sm_next(scatter_cut,:);
    ray_interfaces(num_scatters).ray_index = round(ray_index);
    ray_interfaces(num_scatters).surface_index = round(six_next(scatter_cut));
    ray_interfaces(num_scatters).distance_traveled = l_next(scatter_cut);
    ray_interfaces(num_scatters).n_incident = n_next(scatter_cut,1);
    ray_interfaces(num_scatters).n_transmitted = n_next(scatter_cut,2);
    
    %% get set for next iteration
    % follow reflected and refracted rays that are above the follow_threshold
    refracted_rays_to_follow = refracted_rays(:,7) > follow_threshold(1);
    reflected_rays_to_follow = reflected_rays(:,7) > follow_threshold(2);
    
    p_start = [ p_next(refracted_rays_to_follow,:) ; ...
        p_next(reflected_rays_to_follow,:) ];
    incoming_rays = [ refracted_rays(refracted_rays_to_follow,:) ; ...
        reflected_rays(reflected_rays_to_follow,:) ];

    % identify reflected rays with a negative ray_index (refracted rays
    % also inhered the negative index if they have previously been
    % reflected)
    ray_index = [ ray_index(refracted_rays_to_follow) ; ...
        -abs(ray_index(reflected_rays_to_follow)) ];
    
end
