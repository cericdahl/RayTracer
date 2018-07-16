% function surface_list = CreteMiXGeometry(geospecs)
%
% This is a sample function demonstrating how to build a simple geometry
% for use in RayTracer2 -- this is called by the example script
% RunMiXGeometry
%
% CED 

function surface_list = CreateMiXGeometry(geospecs)

%% fist set default parameters (these can be overwritten by geospecs)
% all length dimenstions in cm

% dimensions
tpc_height = 10;
tpc_rad = 2;
pmtwin_thick = 0.3;
pmt_rad = 1.5*2.54;

% material bulk properties (indices of refraction, scattering lengths)
n_xenon = 1.69;
n_ptfe = 1.3;
n_quartz = 1.59;
n_ss = inf;
rayleigh_xenon = 29;
abslength_xenon = 300;

% surface properties
ptfe_ref = 0.95;
ptfe_siga = 0;
ptfe_Csl = 0;
ptfe_Css = 1; % specular reflection of dielectric interface
ptfe_Cbs = 0;
ptfe_abs = 0; % _abs and _ref are implemented differently -- see documentation

ss_ref = 0; % irrelevant for metallic reflector (n_ss=inf)
ss_siga = 0.2; % surface roughness (rms in radians from normal)
ss_Csl = 1; % reflection off microfact
ss_Css = 0;
ss_Cbs = 0;
ss_abs = .5; % _abs and _ref are implemented differently -- see documentation

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

%% Now would be the place to calculate some derivative quantities useful in geometry building
% not really necessary in a geometry this simple, but will do some anyway
ptfe_unifiedparams = [ptfe_siga, ptfe_ref, ptfe_Csl, ptfe_Css, ptfe_Cbs];

ss_unifiedparams = [ss_siga, ss_ref, ss_Csl, ss_Css, ss_Cbs];

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
    [0, 0, 0], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) < pmt_rad^2, size(p,1), [] ));
surface_list(end).n_outside = n_quartz;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = inf; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'quartz-xenon interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, pmtwin_thick], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) < pmt_rad^2, size(p,1), [] ));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0; % use unified reflectivity instead
surface_list(end).abslength_outside = abslength_xenon; % photons shouldn't get there
surface_list(end).abslength_inside = inf; % xenon
surface_list(end).rayleigh_outside = rayleigh_xenon; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'ss-xenon interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, tpc_height], [0, 0, -1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) < tpc_rad^2, size(p,1), [] ));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_ss;
surface_list(end).surface_type = 'unified';
surface_list(end).absorption = ss_abs; % use unified reflectivity instead
surface_list(end).abslength_outside = abslength_xenon; % photons shouldn't get there
surface_list(end).abslength_inside = inf; % xenon
surface_list(end).rayleigh_outside = rayleigh_xenon; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = ss_unifiedparams;

surface_list(end+1).description = 'ptfe-xenon interface';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], tpc_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) > pmtwin_thick) & (p(:,3,:) < tpc_height), size(p,1), [] ));
surface_list(end).n_outside = n_ptfe;
surface_list(end).n_inside = n_xenon;
surface_list(end).surface_type = 'unified';
surface_list(end).absorption = ptfe_abs; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = abslength_xenon; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = rayleigh_xenon;
surface_list(end).unifiedparams = ptfe_unifiedparams;

surface_list(end+1).description = 'pmtwindow, side wall -- just making this absorbing for now';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], pmt_rad);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) < pmtwin_thick) & (p(:,3,:) > 0), size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_quartz;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1; % use unified reflectivity instead
surface_list(end).abslength_outside = inf; % photons shouldn't get there
surface_list(end).abslength_inside = inf; % xenon
surface_list(end).rayleigh_outside = inf; % photons shouldn't get there
surface_list(end).rayleigh_inside = inf;
surface_list(end).unifiedparams = zeros(1,5);

