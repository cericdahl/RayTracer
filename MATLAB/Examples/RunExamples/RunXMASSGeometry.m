% Sample code to run RayTracer2 on a simple geometry for light collection
% calcs

%% make geospecs structure
geospecs = struct();

geospecs.n_xenon = 1.5426;%1.62;%1.37;%1.62;


%% Create surface list (call CreateMiXGeometry)
surface_list = CreateXMASStempgeometry(geospecs);

%% Create initial set of rays to traces
n_rays = 1e5;

ray_startingpoints = repmat([0, 0, .1], n_rays, 1);
ray_startingpoints(:,3) = -.3*(2.88/2.5826)*log(rand(n_rays,1));
ray_startingpoints(:,1) = .8*sqrt(rand(n_rays,1));
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
