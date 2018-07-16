% function [intersection_points surface_normals distance_traveled crossing_into] = ...
%     RayToCylinder(starting_points, incoming_directions, cylinder_center, cylinder_axis, cylinder_radius)
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
%       cylinder_center         -   3-element vector, giving a point along
%                                     the cylinder axis
%       cylinder_axis           -   3-element vector pointing along the
%                                     cylinder axis
%       cylinder_radius         -   scalar giving cylinder radius
%                                    
%
% 12/15/09, CED

function [intersection_points surface_normals distance_traveled crossing_into] = ...
    RayToCylinder(starting_points, incoming_directions, cylinder_center, cylinder_axis, cylinder_radius)

intersection_points = [];
surface_normals = [];
distance_traveled = [];
crossing_into = [];

%% check inputs
if nargin<4 || numel(cylinder_center)~=3 || numel(cylinder_axis)~=3 || numel(cylinder_radius)~=1 || size(starting_points,2)~=3 || ...
        size(incoming_directions,2)~=3 || size(starting_points,1)~=size(incoming_directions,1)
    disp('Impropper input to RayToCylinder');
    return
end
cylinder_center = cylinder_center(:)';
cylinder_axis = cylinder_axis(:)';
numrays = size(starting_points,1);

%% normalize directions
goodray_cut = sum(incoming_directions.^2,2)>0;
if any(goodray_cut)
    incoming_directions(goodray_cut,:) = incoming_directions(goodray_cut,:) ./ ...
        repmat(abs(sqrt(sum(incoming_directions.^2,2))),1,3);
end

if sum(cylinder_axis.^2) > 0
    cylinder_axis = cylinder_axis ./ abs(sqrt(sum(cylinder_axis.^2)));
else
    disp('Invalid cylinder for RayToCylinder');
    return
end

%% solve quadratic for distance_traveled
% || u + l*v || = r
% a*l^2 + b*l + c = 0, a = 1

x = starting_points - repmat(cylinder_center,numrays,1);
u = (x * cylinder_axis' * cylinder_axis) - x;


v = (incoming_directions * cylinder_axis' * cylinder_axis) - incoming_directions;

a = sum(v .* v, 2);
b = 2 .* sum(u .* v, 2);
c = sum(u .* u, 2) - cylinder_radius.^2;

linear_cut = a==0 & b~=0;
quad_cut = a~=0;
distance_traveled = zeros(numrays,2);
distance_traveled(~(linear_cut | quad_cut),:) = NaN;
if any(linear_cut)
    distance_traveled(linear_cut,:) = ...
        repmat(-c(linear_cut) ./ b(linear_cut), 1, 2);
end
if any(quad_cut)
    distance_traveled(quad_cut,:) = ...
        repmat(-.5 .* b(quad_cut) ./ a(quad_cut),1,2) + ...
        (.5 .* sqrt(b(quad_cut).^2 - 4 .* a(quad_cut) .* c(quad_cut)) ./ a(quad_cut)) * [1 -1];
end

%% find intersection_points
intersection_points = repmat(starting_points,[1,1,2]) + ...
    repmat(reshape(distance_traveled,[],1,2),[1,3,1]) .* repmat(incoming_directions,[1,1,2]);

%% find surface_normals
x = intersection_points(:,:,1) - repmat(cylinder_center,numrays,1);
u = (x * cylinder_axis' * cylinder_axis) - x;
y = intersection_points(:,:,2) - repmat(cylinder_center,numrays,1);
v = (y * cylinder_axis' * cylinder_axis) - y;
surface_normals = zeros(size(intersection_points));
surface_normals(:,:,1) = u ./ cylinder_radius;
surface_normals(:,:,2) = v ./ cylinder_radius;
crossing_into = round(sign(sum(repmat(incoming_directions,[1,1,2]) .* surface_normals,2)));
surface_normals = -surface_normals .* repmat(crossing_into,[1 3 1]);
crossing_into = reshape(crossing_into,[],2);
