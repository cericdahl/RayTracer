% Sample code to run RayTracer2 on a simple geometry for light collection
% calcs

%% make geospecs structure
geospecs = struct();
which_geospecs = 'xebc';

switch which_geospecs
    
    case 'xebc'
        
        %%% jar dimensions
        
        geospecs.itube_ID = 2;
        geospecs.itube_OD = 2.3;
        geospecs.otube_ID = 2.4;
        geospecs.otube_OD = 3;
        geospecs.icap_thick = .5;
        geospecs.ocap_thick = .5;
        
        geospecs.target_height = 1.25*2.54-.5;
end

surface_list = CreateSBCGeometry(geospecs);

%% create interaction list
n_points = 1e4;

zd = rand(n_points, 1) * geospecs.target_height;
r2d = rand(n_points, 1) * .25 * geospecs.otube_ID^2;
rd = sqrt(r2d);
yd = rd;
xd = zeros(size(zd));


%% create spatial bins
pde = zeros(n_points, 1);

for i_p=1:n_points
%     for i_z=1:n_zbins
        
        %% Create initial set of rays to traces
        n_rays = 1e5;
        
        ray_startingpoints = repmat([0, rd(i_p), zd(i_p)], n_rays, 1);
%         ray_startingpoints = [zeros(n_rays, 1), ...
%             sqrt(r2bin_edges(i_r) + diff(r2bin_edges(i_r:(i_r+1)))*rand(n_rays, 1)), ...
%             zbin_edges(i_z) + diff(zbin_edges(i_z:(i_z+1)))*rand(n_rays, 1)];
        
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
        fprintf(1, '%d:  %.1f, %.1f', i_p, rd(i_p), zd(i_p));
        tic;
        [~, absorption_table, raytable] = RayTracer2(ray_startingpoints, rays, surface_list, ...
            max_scatters, 1e-6, 1e-4, -1, 0, 1, 1);
        toc;
        
        %% and analyze result
        
        if any(reshape(isnan(absorption_table),[],1))
            disp('huh, nan''s...');
            absorption_table(isnan(absorption_table))=0;
        end
        
        total_intensity_traced = sum(reshape(absorption_table(:, 1:4, :, :), [], 1));
        total_intensity_remaining = sum(reshape(absorption_table(end, 5, :, :), [], 1));
        
        if abs(total_intensity_traced + total_intensity_remaining - n_rays) > 1
            disp('Accounting problem in RayTracer2, please report bug.');
        end
        
        intensity_detected_by_numscatters = absorption_table(:, 1, end, 2);
        
        total_intensity_detected = sum(intensity_detected_by_numscatters);
        
        total_bulkabsorption = sum(reshape(absorption_table(:,2,:,:),[],1));
        
        pde(i_p) = total_intensity_detected / total_intensity_traced;
%     end
end

%%
save('~cdahl/sbc_uniform_lightcoll.mat', 'pde', 'xd', 'yd', 'zd', 'rd');


%%
n_colors = 1e3;
c_edges = linspace(0, max(pde), n_colors);
c_list = jet(n_colors);
[~, whichcolor] = histc(pde, c_edges);
figure;
clf;
for i_p=1:length(rd)
    plot(rd(i_p)^2, zd(i_p), 'o','color', c_list(whichcolor(i_p), :), 'markerfacecolor', c_list(whichcolor(i_p),:), 'markersize', 6);
    hold on
end


