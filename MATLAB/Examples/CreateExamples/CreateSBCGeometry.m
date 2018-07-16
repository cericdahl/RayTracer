% function surface_list = CreteMiXGeometry(geospecs)
%
% This is a sample function demonstrating how to build a simple geometry
% for use in RayTracer2 -- this is called by the example script
% RunMiXGeometry
%
% CED 

function surface_list = CreateSBCGeometry(geospecs)

%% fist set default parameters (these can be overwritten by geospecs)
% all length dimenstions in cm

%%% jar dimensions

itube_ID = 2;
itube_OD = 2.3;
otube_ID = 2.4;
otube_OD = 3;

icap_thick = .5;
ocap_thick = .5;

PMTgap = 0.01;
PMTwinthick = 0.1;

piezocover_reflectance = 0.85;

target_height = 1.25*2.54-.5;
%%% rad shield

can_obot = -1.125*2.54*4;

% material bulk properties (indices of refraction, scattering lengths)
n_vacuum = 1;
n_xenon = 1.69;
n_quartz = 1.59;
n_jar = n_quartz;
n_target = n_xenon;
rayleigh_xenon = 29;
abslength_xenon = 300;
abslength_quartz = .083;
abslength_silica = 7;
side_absorb = 0;

%% Now load geospecs -- this'll overwrite the above for those fields present in geospecs
% don't really like this way of doing it (eval statements are bad news in
% general) but it gets the job done

if nargin>0 && isstruct(geospecs)
    fn = fieldnames(geospecs);
    for i_f=1:length(fn)
        if ~isempty(geospecs.(fn{i_f}))
            eval([fn{i_f} ' = geospecs.(fn{i_f});']);
        end
    end
end

%% Now build the surface list
surface_list = struct( ...
    'description', {}, ...
    'intersect_function', {}, ...
    'inbounds_function', {}, ...
    'n_outside', {}, ...
    'n_inside', {}, ...
    'surface_type', {}, ...
    'absorption', {}, ...
    'abslength_outside', {}, ...
    'abslength_inside', {}, ...
    'rayleigh_outside', {}, ...
    'rayleigh_inside', {}, ...
    'unifiedparams', {});

%%% quartz bits
surface_list(end+1).description = 'ID of inner tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*itube_ID);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<-icap_thick) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_quartz;
surface_list(end).abslength_inside = inf;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = inf;

surface_list(end+1).description = 'bottom-side of inner tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 (-icap_thick)], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*itube_ID).^2), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_silica;
surface_list(end).abslength_inside = inf;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = inf;

if piezocover_reflectance == 0
    surface_list(end+1).description = 'reflector under of inner tube cap';
    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
        [0 0 (-icap_thick - PMTgap)], [0 0 1]);
    surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*otube_OD).^2), size(p,1), [] ));
    surface_list(end).n_outside = n_vacuum;
    surface_list(end).n_inside = inf;
    surface_list(end).surface_type = 'normal';
    surface_list(end).absorption = 1;
    surface_list(end).abslength_outside = inf;
    surface_list(end).abslength_inside = inf;
    surface_list(end).rayleigh_outside = inf;
    surface_list(end).rayleigh_inside = inf;
else
    surface_list(end+1).description = 'reflector under of inner tube cap';
    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
        [0 0 (-icap_thick - PMTgap)], [0 0 1]);
    surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*itube_ID).^2), size(p,1), [] ));
    surface_list(end).n_outside = n_vacuum;
    surface_list(end).n_inside = inf;
    surface_list(end).surface_type = 'normal';
    surface_list(end).absorption = 1 - piezocover_reflectance;
    surface_list(end).abslength_outside = inf;
    surface_list(end).abslength_inside = inf;
    surface_list(end).rayleigh_outside = inf;
    surface_list(end).rayleigh_inside = inf;
end

surface_list(end+1).description = 'OD of inner tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*itube_OD);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<0) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = side_absorb;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_quartz;
surface_list(end).rayleigh_outside = rayleigh_xenon;
surface_list(end).rayleigh_inside = inf;

surface_list(end+1).description = 'top-side of inner tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 0], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*itube_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_silica;
surface_list(end).rayleigh_outside = rayleigh_xenon;
surface_list(end).rayleigh_inside = inf;

surface_list(end+1).description = 'ID of outer tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*otube_ID);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<target_height) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = side_absorb;
surface_list(end).abslength_outside = abslength_quartz;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = rayleigh_xenon;

surface_list(end+1).description = 'bottom-side of outer tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 target_height], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*otube_ID).^2), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_silica;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = rayleigh_xenon;

surface_list(end+1).description = 'OD of outer tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], .5*otube_OD);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)<(target_height+ocap_thick)) & (p(:,3,:)>(can_obot-2)), size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = inf;
surface_list(end).abslength_inside = abslength_quartz;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = inf;

surface_list(end+1).description = 'top-side of outer tube cap';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 (target_height + ocap_thick)], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*otube_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = inf;
surface_list(end).abslength_inside = abslength_silica;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = inf;

surface_list(end+1).description = 'bottom-side of PMT window';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 (target_height + ocap_thick + PMTgap)], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*otube_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_vacuum;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_silica;
surface_list(end).abslength_inside = inf;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = inf;

surface_list(end+1).description = 'photocathode';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 0 (target_height + ocap_thick + PMTgap + PMTwinthick)], [0 0 1]);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2+p(:,2,:).^2)<=(.5*otube_OD).^2), size(p,1), [] ));
surface_list(end).n_outside = n_vacuum;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = inf;
surface_list(end).abslength_inside = abslength_silica;
surface_list(end).rayleigh_outside = inf;
surface_list(end).rayleigh_inside = inf;
