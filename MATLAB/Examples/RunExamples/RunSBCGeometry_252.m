% Sample code to run RayTracer2 on a simple geometry for light collection
% calcs

%% L_eff function
Lylist = [ ...
    -14.19, 67.86, 81.03 ; ...
    -6.2, 70.15, 81.21 ; ...
    0, 71.41, 81.61 ; ...
    7.4, 73.69, 82.17 ; ...
    20.94, 78.08, 83.9 ; ...
    36.22, 83.38, 86.77 ; ...
    51.34, 89.62, 90.7 ; ...
    64.74, 95.38, 94.79 ; ...
    85.28, 104, 102.6 ; ...
    105.44, 110.89, 111.83 ; ...
    118.42, 116.21, 116.76 ; ...
    135.91, 120.89, 123.94 ; ...
    149.8, 123.67, 129.28 ; ...
    167.6, 125.98, 135.84 ; ...
    185.15, 127.18, 142 ; ...
    ];
Lylist = Lylist ./ repmat([.5*185.15, 119.75, 119.75], size(Lylist, 1), 1);
L_y_Lindhard = @(Er)(10.^interp1(Lylist(:,1), Lylist(:,2), log10(Er), 'linear', 'extrap'));
L_y_Bezrukov = @(Er)(10.^interp1(Lylist(:,1), Lylist(:,3), log10(Er), 'linear', 'extrap'));

%%
inelastic_yield = @(zaid)( ...
    (zaid==54129) .* (39578*.88/13.7) + ...
    (zaid==54131) .* (80185*.88/13.7) + ...
    (zaid==54133) .* (233221*.88/13.7) + ...
    0);

capture_yield = @(zaid)( ...
    (zaid==54128) .* (39578*.88/13.7) + ...
    (zaid==54130) .* (80185*.88/13.7) + ...
    (zaid==54132) .* (233221*.88/13.7) + ...
    0);

%% whichdir
dirlist = { ...
    '/Users/cdahl/data/SBC-MCNP/Run5/Cf252_1e9', ...
    '/home/cdahl/SBC-MCNP/Run5/Cf252_1e9', ...
    };

for i_d=1:length(dirlist)
    recondir = dirlist{i_d};
    if isdir(recondir)
        break
    end
end

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
        geospecs.piezocover_reflectance = 1;
        geospecs.side_absorb = 1;
end

surface_list = CreateSBCGeometry(geospecs);

%% load mcnp dumn1
% rate_122 = .856 * exp(-1.3068);
% rate_136 = .1068 * exp(-1.1676);

% l_122 = 1/2.9003;
% l_136 = 1/2.1792;

mcnp122 = importdata([recondir filesep 'dumn1']);

%% create interaction list
fissionstarts = find(diff([-1;mcnp122(:,1)])>0);
fissionends = find(diff([mcnp122(:,1);inf])>0);
n_points = length(fissionstarts);

%% create spatial bins
n_phc = zeros(n_points, 10);
max_recoil = zeros(n_points, 1);
tic;
for i_p=1:n_points
    if mod(i_p,10)==0
        fprintf(1,'%d of %d, et %.1f minutes\n', i_p, n_points, toc/60);
    end
%     for i_z=1:n_zbins
        
        thisfission = mcnp122(fissionstarts(i_p):fissionends(i_p),:);
        max_recoil(i_p) = max(thisfission(:,7)*1e3);
        n_rays_by_pos = poissrnd(thisfission(:,7)*1e3.*L_y_Bezrukov(thisfission(:,7)*1e3) + ...
            (abs(thisfission(:,4))==1) .* inelastic_yield(thisfission(:,5)) + ...
            (abs(thisfission(:,4))==0) .* capture_yield(thisfission(:,5)) + ...
             0);

        cumpos_ix = [0 ; cumsum(n_rays_by_pos)];
        
        n_rays = cumpos_ix(end);
        ray_startingpoints = zeros(n_rays, 3);
        for i_pos = 1:length(n_rays_by_pos)
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 1) = thisfission(i_pos, 9);
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 2) = thisfission(i_pos, 10);
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 3) = thisfission(i_pos, 11);
        end
        
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
%         fprintf(1, '%d:  %.1f, %.1f', i_p, rd(i_p), zd(i_p));
%         tic;
        [raytracer_output, absorption_table, raytable] = RayTracer2(ray_startingpoints, rays, surface_list, ...
            max_scatters, 1e-6, 1e-4, -1, 1, 1, 1);
%         toc;
        
        %% and analyze result
        
        pcp = zeros(n_rays, 1);
        for i_s=1:length(raytracer_output)
            cut = abs(raytracer_output(i_s).surface_index) == length(surface_list);
            pcp(raytracer_output(i_s).ray_index(cut)) = raytracer_output(i_s).incoming_ray(cut, 7);
        end
        
        n_phc(i_p,:) = sum(repmat(pcp,1,10) > rand(n_rays, 10), 1);
        
%     end
end

%%
save([recondir filesep 'sbc_252_lightcoll_frommcnp_blackwall_shinybottom.mat'], 'n_phc', 'max_recoil');
return

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


