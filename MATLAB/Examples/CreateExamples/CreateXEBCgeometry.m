% function [surface_list rays ray_startingpoints pixels] = CreateXEBCgeometry(geospecs)
%
% 2/5/2015, CED

function [surface_list, rays, ray_startingpoints, pixels] = CreateXEBCgeometry(geospecs)

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

rays = cell(1,2);
ray_startingpoints = cell(1,2);
pixels = cell(1,2);

%% default geospecs
%%% indices of refraction and dimensions used below
n_target = 1.4;
n_jar = 1.458;
n_window = 1.52;
n_air = 1.00;
n_vacuum = 1;
n_mirror = 0;

%%% jar dimensions

itube_ID = 2;
itube_OD = 2.3;
otube_ID = 2.4;
otube_OD = 3;

icap_thick = .5;
ocap_thick = .5;

target_height = 1.25*2.54-.5;

%%% rad shield

ican_ID = 3.1;
ican_OD = 3 + 2.54/4;
ocan_ID = 2.75*2.54;
ocan_OD = 3*2.54;

cangap_bot = 0;
cangap_top = 1.25*2.54-.5;

can_ibot = -1*2.54;
can_obot = -1.125*2.54;
can_itop = 4*2.54;
can_otop = 4.125*2.54;

canwin_bot = 0;
canwin_top = 1.25*2.54-.5;
canwin_width = pi/6;

%%% mirrors
mirror_pitch = pi/4;
mirror_yaw = pi/6;
mirror_height = 1*2.54;
mirror_bot = -2.54;
mirror_top = 3*2.54;
mirror_width = 3*2.54;

%%% window
win_thick = .5;
win_OD = 4;
win_d = 20;
win_h = 8*2.54;
tube_top = 10*2.54;
tube_bot = 4*2.54;

%%% cameras
cam_x = 0;
cam_y = 0;
cam_z = 1;
cam_f = .8;
cam_barreld = 0;
cam_lenstype = 'theta';
cam_sensorsize = [.1 .1];
cam_resolution = [480 640];

cam_pitch = -pi/2;
cam_yaw = 0;
cam_roll = 0;

%% apply geospecs
fn = fieldnames(geospecs);
for n=1:length(fn)
    if ~isempty(geospecs.(fn{n}))
        eval([fn{n} '=geospecs.(fn{n});']);
    end
end

%% derived dimensions
cam_pixelpitch = cam_sensorsize ./ cam_resolution;

%% surface list
%%% quartz bits
surface_list(end+1).description = 'inside of inner tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*itube_ID);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<-icap_thick) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside of inner tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 (-icap_thick)], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*itube_ID).^2), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside of inner tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*itube_OD);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<0) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside of inner tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 0], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*itube_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside of outer tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*otube_ID);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<target_height) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside of outer tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 target_height], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*otube_ID).^2), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside of outer tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*otube_OD);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<(target_height+ocap_thick)) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside of outer tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 (target_height + ocap_thick)], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*otube_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%%% copper bits
surface_list(end+1).description = 'inside of inner can';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*ican_ID);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,3,:)>(cangap_top)) & (p(:,3,:)<(can_otop)) ) | ...
    ((p(:,3,:)>(can_obot)) & (p(:,3,:)<(cangap_bot)))  ...
    ), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'outside of inner can';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*ican_OD);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,3,:)>(cangap_top)) & (p(:,3,:)<(can_itop)) ) | ...
    ((p(:,3,:)>(can_ibot)) & (p(:,3,:)<(cangap_bot)))  ...
    ), size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'inside of outer can';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*ocan_ID);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,3,:)>(can_ibot)) & (p(:,3,:)<(can_itop)) ) & ...
    ((p(:,3,:)<(canwin_bot)) | (p(:,3,:)>(canwin_top)) | ...
    (p(:,2,:)>0) | (abs(p(:,1,:)./p(:,2,:)) > tan(.5*canwin_width)) )  ...
    ), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'outside of outer can';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*ocan_OD);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,3,:)>(can_obot)) & (p(:,3,:)<(can_otop)) ) & ...
    ((p(:,3,:)<(canwin_bot)) | (p(:,3,:)>(canwin_top)) | ...
    (p(:,2,:)>0) | (abs(p(:,1,:)./p(:,2,:)) > tan(.5*canwin_width)) )  ...
    ), size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'outside top of can';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 can_otop], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_OD).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ican_ID).^2) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'inside top of can';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 can_itop], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_ID).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ican_OD).^2) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'inside bottom of can';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 can_ibot], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_ID).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ican_OD).^2) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'outside bottom of can';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 can_obot], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_ID).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ican_ID).^2) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'top of can gap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 cangap_top], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ican_OD).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ican_ID).^2) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'bottom of can gap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 cangap_bot], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ican_OD).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ican_ID).^2) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'top of can window';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 canwin_top], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_OD).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ocan_ID).^2) & ...
    (p(:,2,:)<0) & (abs(p(:,1,:)./p(:,2,:)) <= tan(.5*canwin_width)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = 'bottom of can window';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 canwin_bot], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_OD).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ocan_ID).^2)  & ...
    (p(:,2,:)<0) & (abs(p(:,1,:)./p(:,2,:)) <= tan(.5*canwin_width)) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = '-x side of can gap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0,0,1], [cos(.5*canwin_width), -sin(.5*canwin_width) 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_OD).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ocan_ID).^2)  & ...
    (p(:,3,:) > canwin_bot) & (p(:,3,:) < canwin_top) & ...
    (p(:,2,:) < 0) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

surface_list(end+1).description = '+x side of can gap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0,0,1], [-cos(.5*canwin_width), -sin(.5*canwin_width) 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*ocan_OD).^2) & ...
    ((p(:,1,:).^2+p(:,2,:).^2)>=(.5*ocan_ID).^2)  & ...
    (p(:,3,:) > canwin_bot) & (p(:,3,:) < canwin_top) & ...
    (p(:,2,:) < 0) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%%% viewport
surface_list(end+1).description = 'top of window';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, -win_d win_h], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+(p(:,2,:)+win_d).^2)<=(.5*win_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_air;
surface_list(end).n_inside = n_window;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'bottom of window';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, -win_d (win_h-win_thick)], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+(p(:,2,:)+win_d).^2)<=(.5*win_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_window;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'window_tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, -win_d 0], [0 0 1], .5*win_OD);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)<tube_top) & ...
    (p(:,3,:)>tube_bot) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%%% mirrors
surface_list(end+1).description = '-x mirror';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, -win_d, mirror_height], [sin(mirror_pitch)*sin(mirror_yaw), sin(mirror_pitch)*cos(mirror_yaw), cos(mirror_pitch)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>mirror_bot) & ...
    (p(:,3,:)<mirror_top) & ...
    (p(:,1,:) < 0) & ...
    (abs(p(:,1,:)) < mirror_width) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = n_mirror;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = '-x mirror';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, -win_d, mirror_height], [-sin(mirror_pitch)*sin(mirror_yaw), sin(mirror_pitch)*cos(mirror_yaw), cos(mirror_pitch)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:)>mirror_bot) & ...
    (p(:,3,:)<mirror_top) & ...
    (p(:,1,:) >= 0) & ...
    (abs(p(:,1,:)) < mirror_width) ), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = n_mirror;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%% create ray lists
[raydirections, pixelmap] = GenerateRaysFromCamera(cam_resolution, cam_pixelpitch, .5*(1+cam_resolution), cam_f, ...
    cam_pitch, cam_yaw, cam_roll, cam_barreld, cam_lenstype);

rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
pixels{1} = pixelmap;
ray_startingpoints{1} = repmat([cam_x cam_y-win_d cam_z+win_h],size(raydirections,1),1);
