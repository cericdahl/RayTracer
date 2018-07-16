% function surface_list = CreteMiXGeometry(geospecs)
%
% This is a sample function demonstrating how to build a simple geometry
% for use in RayTracer2 -- this is called by the example script
% RunMiXGeometry
%
% CED 

function surface_list = CreateXMASStempgeometry(geospecs)

%% fist set default parameters (these can be overwritten by geospecs)
% all length dimenstions in cm

% dimensions
tpc_height = 5.5;
pmtwin_thick = 2.6;
pmt_rad = .8;

% material bulk properties (indices of refraction, scattering lengths)
n_xenon = 1.69;
n_mgf2 = 1.44;
abslength_mgf2 = 14.6;

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

surface_list(end+1).description = 'PMT face';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, -pmtwin_thick], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) < pmt_rad^2, size(p,1), [] ));
surface_list(end).n_outside = n_mgf2;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1; % use unified reflectivity instead
surface_list(end).abslength_outside = abslength_mgf2; % photons shouldn't get there
surface_list(end).abslength_inside = inf; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'quartz-xenon interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, 0], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) < pmt_rad^2, size(p,1), [] ));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_mgf2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = abslength_mgf2; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'quartz-ss interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, 0], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) > pmt_rad^2, size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_mgf2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = abslength_mgf2; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'quartz-ss interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, -1], [0, 0, -1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) > pmt_rad^2, size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_mgf2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = abslength_mgf2; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'ss-xenon interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, tpc_height], [0, 0, -1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) < pmt_rad^2, size(p,1), [] ));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = inf; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'ss-xenon interface';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], pmt_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) > 0) & (p(:,3,:) < tpc_height), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_xenon;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = inf; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'pmtwindow, side wall -- just making this absorbing for now';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], pmt_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) < -1) & (p(:,3,:) > -pmtwin_thick), size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_mgf2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = abslength_mgf2; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'pmtwindow, side wall -- just making this absorbing for now';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], pmt_rad+1);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) < 0) & (p(:,3,:) > -1), size(p,1), [] ));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_mgf2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = abslength_mgf2; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);
