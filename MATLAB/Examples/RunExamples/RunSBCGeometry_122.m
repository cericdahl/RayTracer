% Sample code to run RayTracer2 on a simple geometry for light collection
% calcs

%% whichdir
dirlist = { ...
    '/Users/cdahl/data/SBC-MCNP/Run5/Co57_low', ...
    '/home/cdahl/SBC-MCNP/Run5/Co57_low', ...
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
        geospecs.piezocover_reflectance = 0;
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
gammastarts = find(diff([-1;mcnp122(:,1)])>0);
gammaends = find(diff([mcnp122(:,1);inf])>0);
n_points = length(gammastarts);

%% create spatial bins
n_phc = zeros(n_points, 10);
tic;
for i_p=1:n_points
    if mod(i_p,10)==0
        fprintf(1,'%d of %d, et %.1f minutes\n', i_p, n_points, toc/60);
    end
%     for i_z=1:n_zbins
        
        thisgamma = mcnp122(gammastarts(i_p):gammaends(i_p),:);
        n_rays_by_pos = poissrnd(thisgamma(:,7)*1e6*0.88/13.7);
        cumpos_ix = [0 ; cumsum(n_rays_by_pos)];
        
        n_rays = cumpos_ix(end);
        ray_startingpoints = zeros(n_rays, 3);
        for i_pos = 1:length(n_rays_by_pos)
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 1) = thisgamma(i_pos, 9);
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 2) = thisgamma(i_pos, 10);
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 3) = thisgamma(i_pos, 11);
            
            if (thisgamma(i_pos,4) == 0) && (thisgamma(i_pos, 7) > 0.03456)
                xabs_r = -0.0465 * log(rand(1));
                xabs_phi = 2*pi*rand(1);
                xabs_theta = acos(2*rand(1)-1);
                xabs_pos = thisgamma(i_pos, 9:11) + [ ...
                    xabs_r * sin(xabs_theta) * cos(xabs_phi), ...
                    xabs_r * sin(xabs_theta) * sin(xabs_phi), ...
                    xabs_r * cos(xabs_theta) ];
                n_xray = binornd(cumpos_ix(i_pos+1) - cumpos_ix(i_pos), ...
                    0.0298 / thisgamma(i_pos, 7));
                
                ray_startingpoints((1:n_xray) + cumpos_ix(i_pos), :) = ...
                    repmat(xabs_pos, n_xray, 1);
            end
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
save([recondir filesep 'sbc_122_lightcoll_frommcnp_blackwall_xrays.mat'], 'n_phc');
% RunSBCGeometry_252
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


