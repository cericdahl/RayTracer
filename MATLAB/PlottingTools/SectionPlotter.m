function linehandles = SectionPlotter(surface_list,centerpoint,plane_normal,pitchnum,fignum)

if nargin<5 || isempty(fignum)
    figure;
else
    figure(fignum);
end

if nargin<4 || isempty(pitchnum)
    pitchnum = 20;
end

if nargin<3 || isempty(plane_normal)
    plane_normal = [0 0 1];
end

if nargin<2 || isempty(centerpoint)
    centerpoint = [0 0 0];
end

linehandles = zeros(1,length(surface_list));

%% make raylist
starting_points = repmat(centerpoint(:)',pitchnum+1,1);

plane_normal = plane_normal ./ sqrt(sum(plane_normal(:).^2));

if plane_normal(1)==1
    xdir = [0 1 0];
    ydir = [0 0 1];
else
    ydir = cross(plane_normal(:)', [1 0 0]);
    ydir = ydir ./ sqrt(sum(ydir.^2));
    xdir = cross(ydir, plane_normal(:)');
    xdir = xdir ./ sqrt(sum(xdir.^2));
end
   
td = linspace(0,2*pi,pitchnum+1);

rays = cos(td(:))*xdir + sin(td(:))*ydir;

%% make a figure!
for n=1:length(surface_list)
    [p_intersect, s_normal, l_ray, s_orientation] = ...
        surface_list(n).intersect_function(starting_points,rays);
    
    valid_intersection = ...
        ( surface_list(n).inbounds_function(p_intersect) ) & ...
        ( imag(l_ray)==0 ) & ...
        ( s_orientation ~= 0 ) & ...
        ( ~isnan(l_ray) ) & ...
        ( l_ray < inf ) & ...
        ( l_ray > -inf );
    
    inward_points = reshape(valid_intersection & s_orientation>0,[],1);
    outward_points = reshape(valid_intersection & s_orientation<0,[],1);
    
    all_points = reshape(permute(p_intersect,[1 3 2]),[],3);
    xd = sum(all_points .* repmat(xdir,size(all_points,1),1), 2);
    yd = sum(all_points .* repmat(ydir,size(all_points,1),1), 2);
    
    xd_in = xd;
    xd_in(~inward_points) = NaN;
    yd_in = yd;
    yd_in(~inward_points) = NaN;
    
    xd_out = xd;
    xd_out(~outward_points) = NaN;
    yd_out = yd;
    yd_out(~outward_points) = NaN;

    xd_plot = [xd_in ; NaN ; xd_out];
    yd_plot = [yd_in ; NaN ; yd_out];
    
    linehandles(n) = plot(xd_plot, yd_plot);

    hold on;
end


