% function [intersection_points surface_normals distance_traveled crossing_into] = ...
%     RayToTorus(starting_points, incoming_directions, torus_center, torus_axis, torus_r1, torus_r2)
% 
% One of several RayToXXXXX functions used for RayTracing
%
% All RayToXXXXX functions have the following inputs and outputs in common:
%
%   inputs:
%       starting_points         -   N-by-3 array, where N is the number of
%                                     rays being followed.  Gives the
%                                     starting position for each ray
%       incoming_directions     -   N-by-3 array, gives the
%                                     forward-direction of each ray
%   
%   outputs:
%       intersection_points     -   N-by-3-by-M array of intersection
%                                     points, M varies between functions
%       surface_normals         -   N-by-3-by-M array of surface normals.
%                                     Normals point opposite to
%                                     incoming_directions (dot product < 0)
%       distance_traveled       -   N-by-M array of distances (negative for
%                                     hits behind the starting point, NaN,
%                                     +inf, -inf, or imaginary values are
%                                     non-intersections)
%       crossing_into           -   N-by-M array of +1,0,-1, where +1
%                                     indicates the ray crossing into the
%                                     surface, (surface_normal is aligned
%                                     with outward-normal), -1 indicates
%                                     the ray is leaving the surface
%                                     (surface_normal is anti-aligned with
%                                     outward-normal), and 0 indicates a
%                                     glancing blow
%
% RayToCylinder specifics:
% M = 2, rays that miss the cylinder give imaginary outputs for
% distance_traveled
%
% other inputs:
%       torus_center         -   3-element vector, giving the point at
%                                     the torus center
%       torus_axis           -   3-element vector pointing along the
%                                     torus axis
%       torus_r1             -   scalar giving torus large radius
%                                    
%       torus_r2             -   scalar giving torus small radius (height)
%
%
% 9/3/15, CED

function [intersection_points, surface_normals, distance_traveled, crossing_into] = ...
    RayToTorus(starting_points, incoming_directions, torus_center, torus_axis, torus_r1, torus_r2)

intersection_points = [];
surface_normals = [];
distance_traveled = [];
crossing_into = [];

%% check inputs
if nargin<4 || numel(torus_center)~=3 || numel(torus_axis)~=3 || numel(torus_r1)~=1 || numel(torus_r2)~=1 || size(starting_points,2)~=3 || ...
        size(incoming_directions,2)~=3 || size(starting_points,1)~=size(incoming_directions,1)
    disp('Impropper input to RayToTorus');
    return
end
torus_center = torus_center(:)';
torus_axis = torus_axis(:)';
numrays = size(starting_points,1);

%% normalize directions
goodray_cut = sum(incoming_directions.^2,2)>0;
if any(goodray_cut)
    incoming_directions(goodray_cut,:) = incoming_directions(goodray_cut,:) ./ ...
        repmat(abs(sqrt(sum(incoming_directions.^2,2))),1,3);
end

if sum(torus_axis.^2) > 0
    torus_axis = torus_axis ./ abs(sqrt(sum(torus_axis.^2)));
else
    disp('Invalid torus for RayToTorus');
    return
end

%% solve quadratic for distance_traveled
% ( || u + l*v || - r1 )^2 + ( || y + l*w || )^2 = r2^2
% ||u+lv||^2 + r1^2 + ||y+lw||^2 - r2^2 = 2*r1*||u+lv||
% [ uu + 2luv + llvv + r1^2 + yy + 2lyw + llww - r2^2 ] ^2 = 4 * r1^2 * (uu + 2luv + llvv)
% a4*l^4 + a3*l^3 + a2*l^2 + a1*l + a0 = 0

x = starting_points - repmat(torus_center,numrays,1);
y = (x * (torus_axis') * torus_axis);
u = y - x;
w = (incoming_directions * (torus_axis') * torus_axis);
v = w - incoming_directions;
r1 = torus_r1;
r2 = torus_r2;


a4 = (sum(v.*v, 2) + sum(w.*w, 2)).^2;
a3 = 4 * (sum(v.*v, 2)+sum(w.*w, 2)) .* (sum(u.*v, 2)+sum(y.*w, 2));
a2 = 4*(sum(u.*v, 2) + sum(y.*w, 2)).^2 - 4*r1^2*sum(v.*v, 2) + ...
    2*(sum(v.*v, 2) + sum(w.*w, 2)) .* (sum(u.*u, 2) + r1^2 + sum(y.*y, 2) - r2^2);
a1 = 4 .* (sum(u.*v, 2) + sum(y.*w, 2)) .* (sum(u.*u, 2) + r1^2 + sum(y.*y, 2) - r2^2) - ...
    8*r1^2*sum(u.*v, 2);
a0 = (sum(u.*u, 2) + r1^2 + sum(y.*y, 2) - r2^2).^2 - 4*r1^2*sum(u.*u, 2);

a = [a4, a3, a2, a1, a0];

quartic_cut = a(:,1)~=0;
cubic_cut = ~quartic_cut & a(:,2)~=0;
quad_cut = ~(quartic_cut | cubic_cut) & a(:,3)~=0;
linear_cut = ~(quartic_cut | cubic_cut | quad_cut) & a(:,4)~=0;
nan_cut = ~(quartic_cut | cubic_cut | quad_cut | linear_cut);

distance_traveled = zeros(numrays,4);

distance_traveled(nan_cut,:) = NaN;

if any(linear_cut)
    distance_traveled(linear_cut,1:2) = ...
        repmat(-a(linear_cut,5) ./ a(linear_cut,4), 1, 2);
    distance_traveled(linear_cut,3:4) = NaN;
end

if any(quad_cut)
    distance_traveled(quad_cut,1:2) = ...
        repmat(-.5 .* a(quad_cut,4) ./ a(quad_cut,3),1,2) + ...
        (.5 .* sqrt(a(quad_cut,4).^2 - 4 .* a(quad_cut,3) .* a(quad_cut,5)) ./ a(quad_cut,3)) * [1 -1];
    distance_traveled(quad_cut,3:4) = NaN;
end

if any(cubic_cut)
    for ix = ( find(cubic_cut)' )
        distance_traveled(ix,1:3) = roots(a(ix,2:5));
    end
    distance_traveled(cubic_cut,4) = NaN;
end

if any(quartic_cut)
    for ix = ( find(quartic_cut)' )
        distance_traveled(ix,:) = roots(a(ix,:));
    end
end

%% find intersection_points
intersection_points = repmat(starting_points,[1,1,4]) + ...
    repmat(reshape(distance_traveled,[],1,4),[1,3,1]) .* repmat(incoming_directions,[1,1,4]);

%% find surface_normals
surface_normals = zeros(size(intersection_points));
for i_l=1:4
    x = intersection_points(:,:,i_l) - repmat(torus_center,numrays,1);
    y = (x * (torus_axis') * torus_axis);
    u = x - y;
    surface_normals(:,:,i_l) = (1/r2) * ...
        ( u .* (1 - (r1./repmat(sqrt(sum(u.*u,2)),1,3))) + y);
end
    
crossing_into = round(-sign(sum(repmat(incoming_directions,[1,1,4]) .* surface_normals,2)));
surface_normals = surface_normals .* repmat(crossing_into,[1 3 1]);
crossing_into = reshape(crossing_into,[],4);

