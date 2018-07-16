function linehandles = SurfacePlotter(surface_list,boundbox,pitchnum,fignum)

if nargin<4 || isempty(fignum)
    figure;
else
    figure(fignum);
end

if nargin<3 || isempty(pitchnum)
    pitchnum = 20;
end

if nargin<2 || isempty(boundbox)
    boundbox = [-10 10 ; -10 10 ; -10 10];
end

linehandles = zeros(1,length(surface_list));

%% make raylist
xd = boundbox(1,1):(diff(boundbox(1,:))/pitchnum):boundbox(1,2);
yd = boundbox(2,1):(diff(boundbox(2,:))/pitchnum):boundbox(2,2);
zd = boundbox(3,1):(diff(boundbox(3,:))/pitchnum):boundbox(3,2);

x0 = mean(boundbox(1,:));
y0 = mean(boundbox(2,:));
z0 = mean(boundbox(3,:));

z_starting_points_x = reshape(repmat(xd(:),1,numel(yd)),[],1);
z_starting_points_y = reshape(repmat(yd(:)',numel(xd),1),[],1);
z_starting_points_z = zeros(size(z_starting_points_x)) + z0;
z_starting_points = [z_starting_points_x z_starting_points_y z_starting_points_z];
z_rays = repmat([0 0 1],size(z_starting_points,1),1);

y_starting_points_z = reshape(repmat(zd(:),1,numel(xd)),[],1);
y_starting_points_x = reshape(repmat(xd(:)',numel(zd),1),[],1);
y_starting_points_y = zeros(size(y_starting_points_z)) + y0;
y_starting_points = [y_starting_points_x y_starting_points_y y_starting_points_z];
y_rays = repmat([0 1 0],size(y_starting_points,1),1);

x_starting_points_y = reshape(repmat(yd(:),1,numel(zd)),[],1);
x_starting_points_z = reshape(repmat(zd(:)',numel(yd),1),[],1);
x_starting_points_x = zeros(size(x_starting_points_y)) + x0;
x_starting_points = [x_starting_points_x x_starting_points_y x_starting_points_z];
x_rays = repmat([1 0 0],size(x_starting_points,1),1);

starting_points = [x_starting_points ; y_starting_points ; z_starting_points];
rays = [x_rays ; y_rays ; z_rays];

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
    px = squeeze(p_intersect(:,1,:));
    py = squeeze(p_intersect(:,2,:));
    pz = squeeze(p_intersect(:,3,:));
    if any(valid_intersection(:))
        linehandles(n) = plot3(px(valid_intersection), py(valid_intersection), pz(valid_intersection), ...
            'ob', 'markerfacecolor', 'b', 'markersize',4);
    else
        linehandles(n) = plot3(0,0,0, ...
            'ob', 'markerfacecolor', 'b', 'markersize',4);
    end
    hold on;
end


