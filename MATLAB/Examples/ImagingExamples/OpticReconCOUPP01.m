function line_handles = OpticReconCOUPP01(geospecs, imagemask, axes_handles, default_maskmode)

%% set defaults
if nargin<4 || isempty(default_maskmode)
    default_maskmode = 'small';
end

if nargin<3 || isempty(axes_handles)
    figure;
    axes_handles = [0 0];
    axes_handles(1) = subplot(1,2,1);
    axes_handles(2) = subplot(1,2,2);
end

if numel(axes_handles)==1
    axes_handles = [0 0] + axes_handles;
end

if nargin<2 || isempty(imagemask)
    imagemask = cell(1,2);
end

if nargin<1 || isempty(geospecs)
    geospecs = struct();
end

line_handles = zeros(2,3);

%% get 2LGeometry
[surface_list rays ray_startingpoints pixels] = CreateCOUPP01Geometry(geospecs);

%% loop over cameras
for c=1
    %% get resolution of this camera
    resolution = max(pixels{c},[],1);
    
    if prod(resolution)~=size(pixels{c},1)
        disp('resolution mismatch from CreateCOUPP01Geometry output');
        return
    end
    
    pixtest = [repmat((1:resolution(1))',resolution(2),1) ; reshape(repmat(1:resolution(2),resolution(1),1),[],1)];
    if ~all(pixtest(:) == pixels{c}(:))
        disp('unexpected pixels output from CreateCOUPP01Geometry');
        return
    end
    
    %% apply image mask
    if isempty(imagemask{c})
        imagemask{c} = zeros(resolution);
        switch default_maskmode
            case 'none'
                imagemask{c}(:) = -1;
            case 'small'
                switch c
                    case 1
                        imagemask{c}(:) = -1;
                        imagemask{c}(100:400,215:550) = 0;
                    case 2
                        imagemask{c}(:) = -1;
                        imagemask{c}(100:400,215:550) = 0;
                end
            case 'all'
        end
    end
    
    if ~all(size(imagemask{c})==resolution)
        disp('imagemask does not match camera resolution');
        return
    end
    
    pixel_list = find(imagemask{c}==0);
    
    rays_to_trace = rays{c}(pixel_list,:);
    startingpoints_to_trace = ray_startingpoints{c}(pixel_list,:);
    
    %% run RayTracer
    raytracer_output = RayTracer(startingpoints_to_trace, ...
        rays_to_trace, surface_list, 6, 1e-5, [0 100]);
    raytracer_output(7).ray_index = [];
    
    %% build raylists
    all_raylist = raytracer_output(1).ray_index;
        
    hit_quartz_raylist = raytracer_output(3).ray_index( ...
        (raytracer_output(3).ray_index>0) & ...
        ( (raytracer_output(3).surface_index == 3) | ...
        (raytracer_output(3).surface_index == 5) | ...
        (raytracer_output(3).surface_index == 9) ) );
    
    miss_raylist = all_raylist(~ismember(all_raylist, hit_quartz_raylist));
        
    quartz_only_raylist = raytracer_output(4).ray_index( ...
        (raytracer_output(4).ray_index>0) & ...
        ( (raytracer_output(4).surface_index == -3) | ...
        (raytracer_output(4).surface_index == -5) | ...
        (raytracer_output(4).surface_index == -9) ) );
    
    hit_jar_raylist = raytracer_output(4).ray_index( ...
        (raytracer_output(4).ray_index>0) & ...
        ( (raytracer_output(4).surface_index == 1) | ...
        (raytracer_output(4).surface_index == 4) | ...
        (raytracer_output(4).surface_index == 10) | ...
        (raytracer_output(4).surface_index == 11) | ...
        (raytracer_output(4).surface_index == 2) ) );
    
    through_jar_raylist = raytracer_output(5).ray_index( ...
        (raytracer_output(5).ray_index>0) & ...
        ( (raytracer_output(5).surface_index == -1) | ...
        (raytracer_output(5).surface_index == -6) | ...
        (raytracer_output(5).surface_index == -10) | ...
        (raytracer_output(5).surface_index == -11) | ...
        (raytracer_output(5).surface_index == -2) | ...
        (raytracer_output(5).surface_index == -4) ) );
    
    quartz_reflection_raylist = hit_jar_raylist(~ismember(hit_jar_raylist, through_jar_raylist));

    fidmark_raylist = [];
    for n=3:length(raytracer_output)
        fidmark_raylist = [fidmark_raylist ; ...
            raytracer_output(n).ray_index( ...
            (raytracer_output(n).ray_index>0) & ...
            (abs(raytracer_output(n).surface_index)==9) ) ];
    end
    fidmark_raylist = unique(fidmark_raylist);
    
    testmark1_raylist = [];
    for n=4:length(raytracer_output)
        testmark1_raylist = [testmark1_raylist ; ...
            raytracer_output(n).ray_index( ...
            (raytracer_output(n).ray_index>0) & ...
            (abs(raytracer_output(n).surface_index)==10) ) ];
    end
    testmark1_raylist = unique(testmark1_raylist);
    testmark1_raylist = testmark1_raylist(ismember(testmark1_raylist, through_jar_raylist));
    
    testmark2_raylist = [];
    for n=4:length(raytracer_output)
        testmark2_raylist = [testmark2_raylist ; ...
            raytracer_output(n).ray_index( ...
            (raytracer_output(n).ray_index>0) & ...
            (abs(raytracer_output(n).surface_index)==11) ) ];
    end
    testmark2_raylist = unique(testmark2_raylist);
    testmark2_raylist = testmark2_raylist(ismember(testmark2_raylist, through_jar_raylist));

    cf3i_raylist = raytracer_output(5).ray_index( ...
        (raytracer_output(5).ray_index>0) & ...
        ( (raytracer_output(5).surface_index == -1) | ...
        (raytracer_output(5).surface_index == -4) ) );
    cf3i_raylist = cf3i_raylist(~ismember(cf3i_raylist, [testmark1_raylist ; testmark2_raylist]));
    
    h2o_raylist = raytracer_output(5).ray_index( ...
        (raytracer_output(5).ray_index>0) & ...
        (raytracer_output(5).surface_index == -2) );
    h2o_raylist = h2o_raylist(~ismember(h2o_raylist, [testmark1_raylist ; testmark2_raylist]));
    
    
    interface_raylist = raytracer_output(5).ray_index( ...
        (raytracer_output(5).ray_index>0) & ...
        (raytracer_output(5).surface_index == -6) );
    interface_raylist = interface_raylist(~ismember(interface_raylist, [testmark1_raylist ; testmark2_raylist]));
    
    %% fill out imagemask
    imagemask{c}(pixel_list(miss_raylist)) = 1;
    imagemask{c}(pixel_list(quartz_only_raylist)) = 2;
    imagemask{c}(pixel_list(quartz_reflection_raylist)) = 4;
    imagemask{c}(pixel_list(cf3i_raylist)) = 8;
    imagemask{c}(pixel_list(interface_raylist)) = 16;
    imagemask{c}(pixel_list(h2o_raylist)) = 32;
    imagemask{c}(pixel_list(testmark1_raylist)) = 64;
    imagemask{c}(pixel_list(testmark2_raylist)) = 128;
    imagemask{c}(pixel_list(fidmark_raylist)) = 256;

    if any(imagemask{c}(:)==0)
        disp('whoops!  Missed a pixel in the imagemask!');
        return
    end
    
    %% draw boundaries
    vert_boundaries = diff(imagemask{c},1,1)~=0;
    horz_boundaries = diff(imagemask{c},1,2)~=0;
    
    xd = (2:resolution(1)) - .5;
    ydl = (1:resolution(2)) - .5;
    ydu = (1:resolution(2)) + .5;

    yd = (2:resolution(2)) - .5;
    xdl = (1:resolution(1)) - .5;
    xdu = (1:resolution(1)) + .5;
    
    xd = reshape(repmat(xd(:),resolution(2),1), size(vert_boundaries));
    ydl = reshape(repmat(ydl(:)',resolution(1)-1,1), size(vert_boundaries));
    ydu = reshape(repmat(ydu(:)',resolution(1)-1,1), size(vert_boundaries));
    
    yd = reshape(repmat(yd(:)',resolution(1),1), size(horz_boundaries));
    xdl = reshape(repmat(xdl(:),resolution(2)-1,1), size(horz_boundaries));
    xdu = reshape(repmat(xdu(:),resolution(2)-1,1), size(horz_boundaries));
    
    xdata = [ reshape([reshape(xd(vert_boundaries),1,[]) ; reshape(xd(vert_boundaries),1,[]) ; repmat(NaN,1,sum(vert_boundaries(:)))] , [], 1) ;...
        reshape([reshape(xdl(horz_boundaries),1,[]) ; reshape(xdu(horz_boundaries),1,[]) ; repmat(NaN,1,sum(horz_boundaries(:)))], [], 1)];

    ydata = [ reshape([reshape(ydl(vert_boundaries),1,[]) ; reshape(ydu(vert_boundaries),1,[]) ; repmat(NaN,1,sum(vert_boundaries(:)))], [], 1) ; ...
        reshape([reshape(yd(horz_boundaries),1,[]) ; reshape(yd(horz_boundaries),1,[]) ; repmat(NaN,1,sum(horz_boundaries(:)))], [], 1)];
    
    
    for n_a = 1:length(axes_handles)
        axes(axes_handles(n_a));
        linehandle = plot(axes_handles(n_a), xdata, ydata, 'r');
        if ~isempty(linehandle)
            line_handles(n_a,1) = linehandle;
        end
        hold on;
    end

    vert_boundaries = diff(imagemask{c}>0,1,1)~=0;
    horz_boundaries = diff(imagemask{c}>0,1,2)~=0;

    xd = (2:resolution(1)) - .5;
    ydl = (1:resolution(2)) - .5;
    ydu = (1:resolution(2)) + .5;

    yd = (2:resolution(2)) - .5;
    xdl = (1:resolution(1)) - .5;
    xdu = (1:resolution(1)) + .5;

    xd = reshape(repmat(xd(:),resolution(2),1), size(vert_boundaries));
    ydl = reshape(repmat(ydl(:)',resolution(1)-1,1), size(vert_boundaries));
    ydu = reshape(repmat(ydu(:)',resolution(1)-1,1), size(vert_boundaries));

    yd = reshape(repmat(yd(:)',resolution(1),1), size(horz_boundaries));
    xdl = reshape(repmat(xdl(:),resolution(2)-1,1), size(horz_boundaries));
    xdu = reshape(repmat(xdu(:),resolution(2)-1,1), size(horz_boundaries));

    xdata = [ reshape([reshape(xd(vert_boundaries),1,[]) ; reshape(xd(vert_boundaries),1,[]) ; repmat(NaN,1,sum(vert_boundaries(:)))] , [], 1) ;...
        reshape([reshape(xdl(horz_boundaries),1,[]) ; reshape(xdu(horz_boundaries),1,[]) ; repmat(NaN,1,sum(horz_boundaries(:)))], [], 1)];

    ydata = [ reshape([reshape(ydl(vert_boundaries),1,[]) ; reshape(ydu(vert_boundaries),1,[]) ; repmat(NaN,1,sum(vert_boundaries(:)))], [], 1) ; ...
        reshape([reshape(yd(horz_boundaries),1,[]) ; reshape(yd(horz_boundaries),1,[]) ; repmat(NaN,1,sum(horz_boundaries(:)))], [], 1)];

    for n_a = 1:length(axes_handles)
        linehandle = plot(axes_handles(n_a), xdata, ydata, 'b');
        if ~isempty(linehandle)
            line_handles(n_a,2) = linehandle;
        end
    %     line_handles(c,3) = plot(axes_handles(c), test_pixels{c}(:,1), test_pixels{c}(:,2), 'or','markerfacecolor','r','linewidth',1,'markersize',4);

        axis([1 resolution(1) 1 resolution(2)]);
        set(gca,'ydir','reverse');
    end
end


