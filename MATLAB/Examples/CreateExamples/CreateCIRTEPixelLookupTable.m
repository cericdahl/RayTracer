% function Create2LPixelLookupTable(filename)
%
% This function creates a lookup table for the primary ray-trace for each
% pixel in the COUPP 2L chamber.  This lookup table is loaded by
% COUP2L_XYZlookup and passed to FindClosestApproach.
%
% 12/13/2010, CED

function CreateCIRTEPixelLookupTable(filename)

if nargin<1
    filename = 'CIRTE_PixelLookupTable.mat';
end

%% create geospecs
% note that these geospecs have cameras 0 and 1 switched relative to the
% 4kg DAQ -- this is fixed in this function so that the output lookup table
% has camera numbering consistent with the DAQ
geospec_list = { ...
    'n_CF3I', 'n_H2O', 'n_quartz', 'n_glycol', 'n_air', 'n_glass', ...
    'cf3i_mass', 'cf3i_density', ...
    'cam0_focallength', 'cam1_focallength', 'cam0_distortion', 'cam1_distortion', ...
    'cam0_x', 'cam0_y', 'cam0_z', 'cam1_x', 'cam1_y', 'cam1_z', ...
    'cam0_pitch', 'cam0_yaw', 'cam0_roll', 'cam1_pitch', 'cam1_yaw', 'cam1_roll', ...
    'jar_cylrad', 'jar_axrad', 'jar_cylthick', 'jar_axthick', ...
    'jar_pitch', 'jar_yaw', 'jar_roll', ...
    'window_inside', 'window_thickness', ...
    'fid_mark_z1', 'fid_mark_z2', 'fid_mark_rphi', 'fid_mark_length', 'fid_mark_pen', ...
    'surface_test_cyl_z', 'surface_test_cyl_phi', 'surface_test_sph_z', 'surface_test_sph_phi', 'testmark_radius', ...
    };

geospec_defaults = [ ...
    1.31, 1.33, 1.458, 1.33, 1.00, 1.491, ...
    14.7, 2, ...
    .8, .8, 0.0, 0.0, ...
    0, -9.0, 1.65, 0, -9.3, 8.45, ...
    11, 1.1, 0.7, -15, 0.2, 0.7, ...
    .6, .6, .1, .1, ...
    0, 90, -90, ...
    -3.4925, 0.3175, ...
    3.4925, 6.1913, 1.2566, .3, .02, ...
    0, 180, 0, 180, .1];

geospecs = struct();
for n=1:length(geospec_list)
    geospecs.(geospec_list{n}) = geospec_defaults(n);
end

%% Create 2L geometry
[surface_list rays ray_startingpoints pixels] = Create2LGeometry(geospecs);

for c=1:2
    resolution = max(pixels{c},[],1);

    %% run RayTracer
    max_scatters = 10;
    
    raytracer_output = RayTracer(ray_startingpoints{c}, ...
        rays{c}, surface_list, max_scatters, 1e-5, [0 100]);
    
    scatter_points = zeros([prod(resolution), 3, max_scatters]) + inf;
    directions = zeros([prod(resolution), 3, max_scatters]) + inf;
    scatter_here = false(prod(resolution), max_scatters);
    
    for s=1:length(raytracer_output)
        refracted_cut = raytracer_output(s).refracted_ray(:,7) > 0;
        scatter_points(raytracer_output(s).ray_index(refracted_cut),:,s) =  ...
            raytracer_output(s).intersection_point(refracted_cut,:);
        directions(raytracer_output(s).ray_index(refracted_cut),:,s) = ...
            raytracer_output(s).refracted_ray(refracted_cut,1:3);
        scatter_here(raytracer_output(s).ray_index(refracted_cut), s) = true;
        
        reflected_cut = ~refracted_cut & raytracer_output(s).reflected_ray(:,7) > 0;
        scatter_points(raytracer_output(s).ray_index(reflected_cut),:,s) =  ...
            raytracer_output(s).intersection_point(reflected_cut,:);
        directions(raytracer_output(s).ray_index(reflected_cut),:,s) = ...
            raytracer_output(s).reflected_ray(reflected_cut,1:3);
        scatter_here(raytracer_output(s).ray_index(reflected_cut), s) = true;
    end
    
    maxlengths = diff(scatter_points, 1, 3);
    maxlengths = squeeze(sqrt(sum(maxlengths.^2,2)));
    
    num_segments = sum(scatter_here, 2);
    
    pixel_list = struct('starting_point', [], 'direction', [], 'maxlength', []);
    pixel_list(prod(resolution)) = pixel_list;
    
    for p=1:prod(resolution)
        pixel_list(p).starting_point = squeeze(scatter_points(p,:,1:num_segments(p)))';
        pixel_list(p).direction = squeeze(directions(p,:,1:num_segments(p)))';
        pixel_list(p).maxlength = maxlengths(p,1:num_segments(p));
    end
    
    pixel_list = reshape(pixel_list, resolution);
    
    switch c
        case 2
            cam0_pixels = pixel_list;
        case 1
            cam1_pixels = pixel_list;
    end
    
end

save(filename, 'cam0_pixels', 'cam1_pixels', '-mat');

