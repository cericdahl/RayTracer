% function [rays pixels] = GenerateRaysFromCamera(pinhole_position, resolution, ...
%     pixel_pitch, pixel_center, focal_length, pitch, yaw, roll, radial_distortion)
%
% 12/23/09, CED

function [ray_direction pixels] = GenerateRaysFromCamera(resolution, pixel_pitch, ...
    pixel_center, focal_length, pitch, yaw, roll, radial_distortion, lens_type)

%% set defaults
if nargin<9 || isempty(lens_type)
    lens_type = 'tan';
end

if nargin<8 || isempty(radial_distortion)
    radial_distortion = [];
end

if nargin<7 || isempty(roll)
    roll = 0;
end

if nargin<6 || isempty(yaw)
    yaw = 0;
end

if nargin<5 || isempty(pitch)
    pitch = 0;
end

if nargin<4 || numel(focal_length)~=1 || numel(yaw)~=1 || ...
        numel(pitch)~=1 || numel(roll)~=1 || ...
        numel(pixel_center)~=2 || isempty(pixel_pitch) || numel(pixel_pitch)>2 || ...
        numel(resolution)~=2
    disp('impropper input to GenerateRaysFromCamera');
    return
end

if numel(pixel_pitch)==1
    pixel_pitch = [0 0] + pixel_pitch;
end

%% generate rays in camera-frame
% in camera frame, forward is +y, +x is +i, and +z is -j
i_pix = repmat((1:resolution(1))',1,resolution(2));
j_pix = repmat((1:resolution(2)),resolution(1),1);
pixels = [i_pix(:) j_pix(:)];

pixel_location = [(pixel_center(1)-i_pix(:)).*pixel_pitch(1), ...
    zeros(numel(i_pix),1)-focal_length, ...
    -(pixel_center(2)-j_pix(:)).*pixel_pitch(2)];

pixel_d2 = sum(pixel_location(:,[1 3]).^2,2);

effective_f = focal_length * (1 + ...
    sum( repmat(radial_distortion(:)',length(pixel_d2),1) .* ...
    (repmat(focal_length^-2 * pixel_d2,1,length(radial_distortion)).^repmat(1:length(radial_distortion),length(pixel_d2),1)), 2) );

switch lens_type
    case 'theta'
        theta = sqrt(pixel_d2)./effective_f;
    case 'sin'
        theta = asin(sqrt(pixel_d2)./effective_f);
    case 'tan'
        theta = atan(sqrt(pixel_d2)./effective_f);
    otherwise
        theta = atan(sqrt(pixel_d2)./effective_f);
end

phi = atan2(-pixel_location(:,3),-pixel_location(:,1));
ray_direction = [sin(theta).*cos(phi), cos(theta), sin(theta).*sin(phi)];

% ray_direction = ray_direction ./ repmat(abs(sqrt(sum(ray_direction.^2,2))),1,3);

%% rotate camera
M1 = [cos(yaw) -sin(yaw) 0 ; sin(yaw) cos(yaw) 0 ; 0 0 1];
M2 = [1 0 0 ; 0 cos(pitch) -sin(pitch) ; 0 sin(pitch) cos(pitch)];
M3 = [cos(roll) 0 sin(roll) ; 0 1 0 ; -sin(roll) 0 cos(roll)];

M = M1*M2*M3;

ray_direction = (M * (ray_direction'))';

