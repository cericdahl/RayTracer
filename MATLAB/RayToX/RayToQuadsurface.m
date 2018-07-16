% function [intersection_points surface_normals distance_traveled crossing_into] = ...
%     RayToQuadsurface(starting_points, incoming_directions, Q, P, R)
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
% RayToQuadsurface specifics:
% M = 2, rays that miss the sphere give imaginary or NaN outputs for
% distance_traveled
%
% other inputs:
%       Q                       -   3x3 tensor
%       P                       -   3 element vector
%       R                       -   scalar
%
% The surface is that which satisfies the equation
%      x' * Q * x   +   P' * x   +   R   =   0
%
% (where here x is 3x1, Q is 3x3, and P is 3x1)
%
% the outward-pointing normal at x is taken in the direction of
%
%      2 * Q * x   +   P
%
%
% 12/15/09, CED

function [intersection_points surface_normals distance_traveled crossing_into] = ...
    RayToQuadsurface(starting_points, incoming_directions, Q, P, R)

intersection_points = [];
surface_normals = [];
distance_traveled = [];
crossing_into = [];

%% check inputs
if nargin<4 || numel(Q)~=9 || numel(P)~=3 || numel(R)~=1 || size(starting_points,2)~=3 || ...
        size(incoming_directions,2)~=3 || numel(starting_points)~=numel(incoming_directions) || ...
        length(size(starting_points))~=2 || length(size(incoming_directions))~=2
    disp('Impropper input to RayToQuadsurface');
    return
end
Q = reshape(Q,[3 3]);
P = P(:);
numrays = size(starting_points,1);

%% normalize directions
goodray_cut = sum(incoming_directions.^2,2)>0;
if any(goodray_cut)
    incoming_directions(goodray_cut,:) = incoming_directions(goodray_cut,:) ./ ...
        repmat(abs(sqrt(sum(incoming_directions.^2,2))),1,3);
end

%% solve quadratic for distance_traveled
% a*l^2 + b*l + c = 0
a = sum((incoming_directions * Q) .* incoming_directions, 2);
b = incoming_directions * P + ...
    sum((incoming_directions * Q) .* starting_points, 2) + ...
    sum((starting_points * Q) .* incoming_directions, 2);
c = R + (starting_points * P) + ...
    sum((starting_points * Q) .* starting_points, 2);

% linear_cut = a==0 & b~=0;
% the previous line is strictly true, but we also want to avoid rounding error here...
linear_cut = b~=0;
linear_cut(linear_cut) = abs(4.*a(linear_cut).*c(linear_cut)./(b(linear_cut).*b(linear_cut)))<(100*eps);
quad_cut = a~=0 & ~linear_cut;
distance_traveled = zeros(numrays,2);
distance_traveled(~(linear_cut | quad_cut),:) = NaN;
if any(linear_cut)
    distance_traveled(linear_cut,1) = ...
        -c(linear_cut) ./ b(linear_cut);
    distance_traveled(linear_cut,2) = ...
        -b(linear_cut) ./ a(linear_cut);
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
surface_normals = zeros(size(intersection_points));
for n=1:size(surface_normals,3)
    surface_normals(:,:,n) = 2 * intersection_points(:,:,n) * Q + repmat(P', numrays, 1);
    goodnormal_cut = sum(surface_normals(:,:,n).^2, 2)>0;
    surface_normals(goodnormal_cut,:,n) = surface_normals(goodnormal_cut,:,n) ./ ...
        repmat(abs(sqrt(sum(surface_normals(goodnormal_cut,:,n).^2, 2))),1,3);
end
crossing_into = round(-sign(sum(repmat(incoming_directions,[1,1,2]) .* surface_normals,2)));
surface_normals = surface_normals .* repmat(crossing_into,[1 3 1]);
crossing_into = reshape(crossing_into,[],2);
