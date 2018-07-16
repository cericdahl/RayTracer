% function [nearest_approach, D, half_d2D2dl2] = NearestApproach(a, b, c, d)
%
% Finds the point of nearest approach between lines ab and cd
%
% inputs
%           a, b, c, d  -   Nxd arrays, where N is the number of problems
%                             to solve (pairs of lines), and d is the
%                             dimension of the space.
% outputs
%           nearest_approach - Nxd array, giving the nearest approach
%                                between the lines connected a-to-b and
%                                c-to-do
%           D                - Distance of closest approach (full distance
%                                between lines, or twice the distance from
%                                the nearest_approach to either line).  Nx1
%                                vector
%           half_d2D2dl2     - 2nd derivative of the D^2 with respect to
%                                distance traveled along either line, times
%                                .5.  This varies from 0 (parallel lines)
%                                to 1 (perpendicular lines).  Nx1 vector
%
% 12/15/09, CED

function [nearest_approach, D, half_d2D2dl2] = NearestApproach(a, b, c, d)

nearest_approach = [];
D = [];
d2Ddl2 = [];

%% check inputs
if nargin<4 || length(size(a))~=2 || ...
        length(size(a))~=length(size(b)) || ~all(size(a)==size(b)) || ...
        length(size(a))~=length(size(c)) || ~all(size(a)==size(c)) || ...
        length(size(a))~=length(size(d)) || ~all(size(a)==size(d))
    disp('Impropper input to Nearest Approach');
    return
end

numdims = size(a,2);

%% find nearest approaches
u = (c-b) - repmat(sum((c-b).*(a-b),2) ./ sum((a-b).*(a-b),2), 1, numdims) .* (a-b);
v = (d-c) - repmat(sum((d-c).*(a-b),2) ./ sum((a-b).*(a-b),2), 1, numdims) .* (a-b);

length1 = -sum(u.*v,2)./sum(v.*v,2);
nearest_approach_line1 = c + repmat(length1,1,numdims).*(d-c);

u = (a-d) - repmat(sum((a-d).*(c-d),2) ./ sum((c-d).*(c-d),2), 1, numdims) .* (c-d);
v = (b-a) - repmat(sum((b-a).*(c-d),2) ./ sum((c-d).*(c-d),2), 1, numdims) .* (c-d);

length2 = -sum(u.*v,2)./sum(v.*v,2);
nearest_approach_line2 = a + repmat(length2,1,numdims).*(b-a);

nearest_approach = .5 .* (nearest_approach_line1 + nearest_approach_line2);
D = abs(sqrt(sum((nearest_approach_line1-nearest_approach_line2).^2,2)));
half_d2D2dl2 = sum(v.*v,2) ./ sum((a-b).*(a-b),2);
