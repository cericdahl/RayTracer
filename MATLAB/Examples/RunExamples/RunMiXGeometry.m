% Sample code to run RayTracer2 on a simple geometry for light collection
% calcs

%% make geospecs structure
geospecs = struct();

% all the fields in geospecs will overwrite defaults in CreateMiXGeometry.m
% -- called in next cell

% Note, I've guessed at typical values for a lot of this, just meant to be
% an example of how this works

geospecs.tpc_height = 15; % cm
geospecs.n_xenon = 1.69;

ptfe_reflectivity_mode = 'dielectric_plus_diffuse';

switch ptfe_reflectivity_mode
    case 'dielectricspecular_plus_diffuse'  % same as G4's polished-back-painted dielectric-dielectric unified model, and virtually the same as Coimbra's model
        
        geospecs.n_ptfe = 1.3;
        geospecs.ptfe_ref = 0.95; % fraction of transmitted rays that get reflected inside ptfe (always diffuse)
        geospecs.ptfe_siga = 0; % surface roughness (rms in radians from average normal)
        geospecs.ptfe_Csl = 0; % prob of specular reflection off microfacet
        geospecs.ptfe_Css = 1; % prob of specular reflection off average normal
        geospecs.ptfe_Cbs = 0; % prob of retroreflection
        % Csl, Css, and Cbs refer to reflections of the dielectric
        % interface (e.g. they have zero effect if n_ptfe = n_xenon).
        % There is also diffuse reflection off that interface, with
        % probability 1-(Csl+Css+Cbs)
        geospecs.ptfe_abs = 0; % probability of absorption, applied after the case (hits both rays reflected off the dielectric interface
        % and rays undergoing diffuse reflection inside the ptfe
        
    case 'pure_diffuse'
        geospecs.n_ptfe = inf; % all rays reflected off dielectric interface
        geospecs.ptfe_ref = 0; % this parameter now does nothing, since all rays are reflected
        geospecs.ptfe_siga = 0; % this parameter also does nothing in this case, since Csl=0
        geospecs.ptfe_Csl = 0; % prob of specular reflection off microfacet
        geospecs.ptfe_Css = 0; % prob of specular reflection off average normal
        geospecs.ptfe_Cbs = 0; % prob of retroreflection
        % so there's 100% diffuse reflection off the dielectric interface
        geospecs.ptfe_abs = 0.05; % absorption probability
        
    case 'pure_diffuse_again' % should give identical results to previous case (except for polarization effects
        % -- diffuse scattering inside ptfe erases polarization, while
        % diffuse scattering implemented at the dielectric interface treats
        % polarization as if it were a single reflection.
        geospecs.n_ptfe = geospecs.n_xenon; % same as xenon
        geospecs.ptfe_ref = 0.95;
        geospecs.ptfe_abs = 0;
        % and none of the Csl, etc, matter, since all rays are transmitted
        % through the dielectric interface
end

% Lots of other options you can set, see CreateMiXGeometry.m for rest of
% geometry, and expand as needed.


%% Create surface list (call CreateMiXGeometry)
surface_list = CreateMiXGeometry(geospecs);

%% Create initial set of rays to traces
n_rays = 1e4;

ray_startingpoints = repmat([0, 0, geospecs.tpc_height-.01], n_rays, 1);

rays = zeros(n_rays, 10);

% each ray starts unpolarized with intensity =1
rays(:, 7) = 1;

% set ray directions (random into 4*pi)
costheta = 1 - 2*rand(n_rays, 1);
sintheta = sqrt(1-costheta.^2);
phi = 2*pi*rand(n_rays,1);
rays(:, 3) = costheta;
rays(:, 1) = sintheta .* cos(phi);
rays(:, 2) = sintheta .* sin(phi);

% set ray polarization reference axis (anything perpendicular to the
% direction, since these are unpolarized)
rays(:, 4:6) = cross(repmat([1, 0, 0], n_rays, 1), rays(:, 1:3));
bad_polref = sum(rays(:, 4:6).^2, 2) == 0;
rays(bad_polref, 4:6) = cross(repmat([0, 1, 0], sum(bad_polref), 1), rays(bad_polref, 1:3));
rays(:, 4:6) = rays(:, 4:6) ./ repmat(abs(sqrt(sum(rays(:, 4:6).^2, 2))), 1, 3);

%% Now run RayTracer2 (here just getting absorption table output)
% Full photon history is in first output -- set last argument to 1 to fill
% out full photon history
max_scatters = 100;
tic;
[~, absorption_table, raytable] = RayTracer2(ray_startingpoints, rays, surface_list, ...
    max_scatters, 1e-4, 0, -1, 0, 1, 1);
toc;

%% and analyze result

total_intensity_traced = sum(reshape(absorption_table(:, 1:4, :, :), [], 1));
total_intensity_remaining = sum(reshape(absorption_table(end, 5, :, :), [], 1));

if abs(total_intensity_traced + total_intensity_remaining - n_rays) > 1e-3
    disp('Accounting problem in RayTracer2, please report bug.');
end

intensity_detected_by_numscatters = absorption_table(:, 1, 1, 1);

total_intensity_detected = sum(intensity_detected_by_numscatters);

total_bulkabsorption = sum(reshape(absorption_table(:,2,:,:),[],1));

total_ptfeabsorption = sum(reshape(absorption_table(:,1,4,:),[],1));

total_ssabsorption = sum(reshape(absorption_table(:,1,3,:),[],1));

pde = total_intensity_detected / total_intensity_traced;
