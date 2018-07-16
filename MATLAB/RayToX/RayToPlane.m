% function [intersection_points surface_normals distance_traveled crossing_into] = ...
%     RayToPlane(starting_points, incoming_directions, plane_point, plane_normal)
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
% M = 1, rays that miss the plane give NaN or Inf outputs for
% distance_traveled, plane_normal defines outward-normal
%
% other inputs:
%       plane_point             -   3-element vector, giving a point on
%                                     the plane
%       plane_normal            -   3-element vector normal to the plane,
%                                     pointing outward
%                                    
%
% 12/15/09, CED

function [intersection_points surface_normals distance_traveled crossing_into] = ...
    RayToPlane(starting_points, incoming_directions, plane_point, plane_normal)

intersection_points = [];
surface_normals = [];
distance_traveled = [];
crossing_into = [];

%% check inputs
if nargin<4 || numel(plane_point)~=3 || numel(plane_normal)~=3 || size(starting_points,2)~=3 || ...
        size(incoming_directions,2)~=3 || size(starting_points,1)~=size(incoming_directions,1)
    disp('Impropper input to RayToPlane');
    return
end
plane_point = plane_point(:)';
plane_normal = plane_normal(:)';
numrays = size(starting_points,1);

%% normalize directions
goodray_cut = sum(incoming_directions.^2,2)>0;
if any(goodray_cut)
    incoming_directions(goodray_cut,:) = incoming_directions(goodray_cut,:) ./ ...
        repmat(abs(sqrt(sum(incoming_directions.^2,2))),1,3);
end

if sum(plane_normal.^2) > 0
    plane_normal = plane_normal ./ abs(sqrt(sum(plane_normal.^2)));
else
    disp('Invalid plane for RayToPlane');
    return
end

%% find intersection points
oldstate = warning('query','MATLAB:divideByZero');
warning('off','MATLAB:divideByZero');
distance_traveled = ((repmat(plane_point,numrays,1) - starting_points) * plane_normal') ./ ...
    (incoming_directions * plane_normal');
warning(oldstate.state,'MATLAB:divideByZero');

intersection_points = starting_points + repmat(distance_traveled,1,3) .* incoming_directions;
surface_normals = -repmat(plane_normal,numrays,1) .* ...
    repmat(sign(incoming_directions * plane_normal'),1,3);
crossing_into = round(-sign(incoming_directions * plane_normal'));
