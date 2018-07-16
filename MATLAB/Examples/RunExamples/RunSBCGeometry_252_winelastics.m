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

L_y_gamma = 1e3*0.88/13.7;

XenonInelasticTableMaker;
% the above creates the following variables:
% xenon_ncapture_q_keV
% xenon_gamma_interactions
% xenon_nnxs_tables
% xenongammatables

%%
% inelastic_yield = @(zaid)( ...
%     (zaid==54129) .* (39578*.88/13.7) + ...
%     (zaid==54131) .* (80185*.88/13.7) + ...
%     (zaid==54133) .* (233221*.88/13.7) + ...
%     0);
% 
% capture_yield = @(zaid)( ...
%     (zaid==54128) .* (39578*.88/13.7) + ...
%     (zaid==54130) .* (80185*.88/13.7) + ...
%     (zaid==54132) .* (233221*.88/13.7) + ...
%     0);

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
        geospecs.piezocover_reflectance = 0;
        geospecs.side_absorb = 0;
end

surface_list = CreateSBCGeometry(geospecs);

%% load mcnp dumn1
% rate_122 = .856 * exp(-1.3068);
% rate_136 = .1068 * exp(-1.1676);

% l_122 = 1/2.9003;
% l_136 = 1/2.1792;

mcnp252 = importdata([recondir filesep 'dumn1']);

%% create interaction list
fissionstarts = find(diff([-1;mcnp252(:,1)])>0);
fissionends = find(diff([mcnp252(:,1);inf])>0);

if 0
    mcnp252 = mcnp252(1:fissionends(1e3),:);
    fissionstarts = find(diff([-1;mcnp252(:,1)])>0);
    fissionends = find(diff([mcnp252(:,1);inf])>0);
end

n_points = length(fissionstarts);
n_lines = fissionends - fissionstarts + 1;

%%
g_cut = mcnp252(:,3)==2;
el_cut = mcnp252(:,3)==1 & mcnp252(:,4)==-99;
nonel_cut = mcnp252(:,3)==1 & mcnp252(:,4)~=-99;

cum_nonel_cut = cumsum(nonel_cut);
no_nonel = diff([0;cum_nonel_cut(fissionends)])==0;

cum_g_cut = cumsum(g_cut);
no_g = diff([0;cum_g_cut(fissionends)])==0;

g_energy = mcnp252(:,7)*1e3;
g_energy(~g_cut) = 0;
cum_g_energy = cumsum(g_energy);
g_energy = diff([0;cum_g_energy(fissionends)]);

nr_energy = mcnp252(:,7)*1e3;
nr_energy(g_cut) = 0;
sorted_nr_energy = sortrows([mcnp252(:,1),nr_energy],[1,2]);
max_nr_energy = sorted_nr_energy(fissionends,2);
next_to_max_nr_energy = zeros(size(max_nr_energy));
next_to_max_nr_energy(n_lines>1) = sorted_nr_energy(fissionends(n_lines>1)-1,2);


%% create spatial bins
n_phc = zeros(n_points, 10);
max_recoil = max_nr_energy;
nexttomax_recoil = next_to_max_nr_energy;
er_energy = g_energy;
elastic_only_cut = no_nonel;
tic;
for i_p=1:n_points
    if mod(i_p,10)==0
        fprintf(1,'%d of %d, et %.1f minutes\n', i_p, n_points, toc/60);
    end
    
    if max_recoil(i_p) < .1 % won't make a bubble so not worth simulating
        continue
    end
    
%     for i_z=1:n_zbins
        
        thisfission = mcnp252(fissionstarts(i_p):fissionends(i_p),:);
%         max_recoil(i_p) = max(thisfission(:,7)*1e3);
        this_gcut = thisfission(:,3)==2;
        n_rays_by_pos = poissrnd( ...
            thisfission(:,7)*1e3.*L_y_Lindhard(thisfission(:,7)*1e3).*(~this_gcut) + ...
            thisfission(:,7)*1e3*L_y_gamma.*(this_gcut) + ...
            0);

        cumpos_ix = [0 ; cumsum(n_rays_by_pos)];
        
        n_rays = cumpos_ix(end);
        ray_startingpoints = zeros(n_rays, 3);
        for i_pos = 1:length(n_rays_by_pos)
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 1) = thisfission(i_pos, 9);
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 2) = thisfission(i_pos, 10);
            ray_startingpoints((cumpos_ix(i_pos)+1):cumpos_ix(i_pos+1), 3) = thisfission(i_pos, 11);
        end
        
        if (~no_nonel(i_p)) && no_g(i_p)
            % now we have to do some work to get some gamma-ish stuff in
            % here
            for i_s=1:size(thisfission,1)
                if thisfission(i_s,3)==2
                    continue
                end
                if thisfission(i_s,4)==-99
                    continue
                end
                if thisfission(i_s,4)==0
                    fn = sprintf('zaid%d',thisfission(i_s,5)+1);
                    if ~isfield(xenon_nnxs_tables, fn)
                        continue
                    else
                        thiscap_ix = find(xenon_ncapture_q_keV(:,1)==thisfission(i_s,5),1,'first');
                        if isempty(thiscap_ix)
                            disp('WHIIPS!');
                            continue
                        end
                    end
                    En = thisfission(i_s,16)*1e3 + xenon_ncapture_q_keV(thiscap_ix,2);
                else
                    fn = sprintf('zaid%d',thisfission(i_s,5));
                    En = thisfission(i_s,16)*1e3;
                end
                % now treat this as an inelastic scatter on zaid given by
                % fn and incoming neutron of energy En
                
                n_exc_states = length(xenon_nnxs_tables.(fn));
                p_exc_states = zeros(1,n_exc_states);
                for i_st=1:n_exc_states
                    p_exc_states(i_st) = interp1( ...
                        [0;log(xenon_nnxs_tables.(fn){i_st}(:,1));log(1e8)], ...
                        [0;xenon_nnxs_tables.(fn){i_st}(:,2);xenon_nnxs_tables.(fn){i_st}(end,2)], ...
                        log(En*1e3), 'linear');
                end
                cdf_exc_states = cumsum(p_exc_states);
                cdf_exc_states = cdf_exc_states / cdf_exc_states(end);
                whichstate = find(cdf_exc_states>rand(1), 1, 'first');
                
                % now get the list of gammas and/or conversion electrons
                IC_energy = 0;
                gammalist = [];
                while whichstate > 0
                    thisroll = rand(1,2);
                    cdf_nextstates = cumsum(xenongammatables.(fn){whichstate,2}(:,2));
                    cdf_nextstates = cdf_nextstates / cdf_nextstates(end);
                    nextstate_ix = find(cdf_nextstates>thisroll(1), 1, 'first');
                    nextstate = xenongammatables.(fn){whichstate,2}(nextstate_ix,1);
                    ICcoeff = xenongammatables.(fn){whichstate,2}(nextstate_ix,3);
                    isIC = (ICcoeff/(1+ICcoeff)) > thisroll(2);
                    this_gamma = xenongammatables.(fn){whichstate,1};
                    whichstate = nextstate;
                    if nextstate > 0
                        this_gamma = this_gamma - xenongammatables.(fn){nextstate,1};
                    elseif nextstate < 0 % indicates long half-life to ground state
                        this_gamma = 0;
                        isIC = true;
                    end
                    if isIC
                        IC_energy = IC_energy + this_gamma;
                    else
                        gammalist(end+1) = this_gamma;
                    end
                end
                   
                gammalist = gammalist(:);
                
                % now see which gammas we get to keep
                if ~isempty(gammalist)
                    gamma_int_dist = 1./(2.5826*exp(interp1(log(xenon_gamma_interactions(:,1)*1e3), log(xenon_gamma_interactions(:,end)), log(gammalist), 'linear', 'extrap')));
                    gamma_int_dist = -gamma_int_dist .* log(rand(size(gamma_int_dist)));
                    gamma_pe_prob = ...
                        exp(interp1(log(xenon_gamma_interactions(:,1)*1e3), log(xenon_gamma_interactions(:,4)), log(gammalist), 'linear', 'extrap') - ...
                        interp1(log(xenon_gamma_interactions(:,1)*1e3), log(xenon_gamma_interactions(:,end)), log(gammalist), 'linear', 'extrap'));
                    is_pe = gamma_pe_prob > rand(size(gamma_pe_prob));
                    costheta = 1-2*rand(length(gammalist),1);
                    phi = 2*pi*rand(length(gammalist),1);
                    gamma_int_disp = repmat(gamma_int_dist(:),1,3) .* [sqrt(1-costheta.^2).*cos(phi), sqrt(1-costheta.^2).*sin(phi), costheta];
                    gamma_int_pos = gamma_int_disp + repmat(thisfission(i_s, 9:11),length(gammalist),1);
                    in_xe = (gamma_int_pos(:,3)) > 0 & (gamma_int_pos(:,3) < geospecs.target_height) & ...
                        (sum(gamma_int_pos(:,1:2).^2,2) < (.5*geospecs.otube_ID));
                else
                    in_xe = false(0,1);
                    gamma_int_pos = zeros(0,3);
                    is_pe = false(0,1);
                end
                
                if IC_energy==0 && ~any(in_xe)
                    continue
                end
                
                % ok, now make photons to follow
                er_pos = [thisfission(i_s, 9:11) ; gamma_int_pos];
                er_kevee = [IC_energy ; gammalist];
                compton_cut = [false ; ~is_pe];
                if any(compton_cut)
                    er_kevee(compton_cut) = er_kevee(compton_cut) .* rand(sum(compton_cut),1);
                end
                in_xe = [(IC_energy>0);in_xe];
                
                er_pos = er_pos(in_xe,:);
                er_kevee = er_kevee(in_xe,:);
                
                
                er_n_rays_by_pos = poissrnd(er_kevee*L_y_gamma);
                er_cumpos_ix = [0 ; cumsum(er_n_rays_by_pos)];
                
                er_n_rays = er_cumpos_ix(end);
                er_ray_startingpoints = zeros(er_n_rays, 3);
                for i_erpos = 1:length(er_n_rays_by_pos)
                    er_ray_startingpoints((er_cumpos_ix(i_erpos)+1):er_cumpos_ix(i_erpos+1), 1) = er_pos(i_erpos, 1);
                    er_ray_startingpoints((er_cumpos_ix(i_erpos)+1):er_cumpos_ix(i_erpos+1), 2) = er_pos(i_erpos, 2);
                    er_ray_startingpoints((er_cumpos_ix(i_erpos)+1):er_cumpos_ix(i_erpos+1), 3) = er_pos(i_erpos, 3);
                end
                
                ray_startingpoints = [ray_startingpoints ; er_ray_startingpoints];
                n_rays = n_rays + er_n_rays;
                er_energy(i_p) = er_energy(i_p) + sum(er_kevee);
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
%         nominal_rays_to_simulate = size(rays,1);
%         max_rays_to_simulate = 1e4;
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
save([recondir filesep 'sbc_252_lightcoll_frommcnp_win_wolonglived.mat'], 'n_phc', 'max_recoil', 'nexttomax_recoil', 'er_energy', 'elastic_only_cut');
return

%%
% n_colors = 1e3;
% c_edges = linspace(0, max(pde), n_colors);
% c_list = jet(n_colors);
% [~, whichcolor] = histc(pde, c_edges);
% figure;
% clf;
% for i_p=1:length(rd)
%     plot(rd(i_p)^2, zd(i_p), 'o','color', c_list(whichcolor(i_p), :), 'markerfacecolor', c_list(whichcolor(i_p),:), 'markersize', 6);
%     hold on
% end


