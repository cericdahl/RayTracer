% function [intersection_points surface_normals distance_traveled crossing_into] = ...
%     RayToSphere(starting_points, incoming_directions, sphere_center, sphere_radius)
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
% RayToSphere specifics:
% M = 2, rays that miss the sphere give imaginary outputs for
% distance_traveled
%
% other inputs:
%       sphere_center           -   3-element vector, giving the center of
%                                     the sphere
%       sphere_radius           -   scalar giving sphere radius
%                                    
%
% 12/15/09, CED

function [intersection_points surface_normals distance_traveled crossing_into] = ...
    RayToSphere(starting_points, incoming_directions, sphere_center, sphere_radius)

intersection_points = [];
surface_normals = [];
distance_traveled = [];
crossing_into = [];

%% check inputs
if nargin<4 || numel(sphere_center)~=3 || numel(sphere_radius)~=1 || size(starting_points,2)~=3 || ...
        size(incoming_directions,2)~=3 || size(starting_points,1)~=size(incoming_directions,1)
    disp('Impropper input to RayToPlane');
    return
end
sphere_center = sphere_center(:)';
numrays = size(starting_points,1);

%% normalize directions
goodray_cut = sum(incoming_directions.^2,2)>0;
if any(goodray_cut)
    incoming_directions(goodray_cut,:) = incoming_directions(goodray_cut,:) ./ ...
        repmat(abs(sqrt(sum(incoming_directions.^2,2))),1,3);
end

%% solve quadratic for distance_traveled
% a*l^2 + b*l + c = 0, a = 1
a = sum(incoming_directions.^2, 2);
b = 2 .* sum(incoming_directions .* (starting_points - repmat(sphere_center,numrays,1)), 2);
c = sum((starting_points - repmat(sphere_center,numrays,1)).^2, 2) - sphere_radius.^2;

distance_traveled = repmat(-.5 .* b ./ a,1,2) + ...
    (.5 .* sqrt(b.^2 - 4 .* a .* c) ./ a) * [1 -1];

%% find intersection_points
intersection_points = repmat(starting_points,[1,1,2]) + ...
    repmat(reshape(distance_traveled,[],1,2),[1,3,1]) .* repmat(incoming_directions,[1,1,2]);

%% find surface_normals
surface_normals = (intersection_points - repmat(sphere_center,[numrays,1,2])) ./ sphere_radius;
crossing_into = round(-sign(sum(repmat(incoming_directions,[1,1,2]) .* surface_normals,2)));
surface_normals = surface_normals .* repmat(crossing_into,[1 3 1]);
crossing_into = reshape(crossing_into,[],2);
