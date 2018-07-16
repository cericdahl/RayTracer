


% all_params = [1 5 1823];
% which_fit_to_fit = logical([1 0 0]);
% 
% fitfun = @(params)OpticsChisqCalc(params,all_params,which_fit_to_fit, pixel_list, pixel_targets, pixel_target_surface, pixel_target_errors);
%
% [best_fit min_chisq] = fminunc(fitfun, all_params(which_fit_to_fit) )


function chisq = OpticsChisqCalc(fitparams, which_fit, rays_list, target)
%    rays_list, target_surfaces, pixel_target_errors
%target with fields: target.(position,surface,error)   
%position currently hard coded.
%rays_list=[126,540;109,887;697,254;693,429;690,629;685,901;910,935]


% %% set default geospecs for COUPP60
%         toplevel = 0;
%         geospecs.n_target = 1.197;
% %         geospecs.n_buffer = 1.33;
% %         geospecs.n_jar = 1.458;
% %         geospecs.n_hydraulic = 1.434;
%         geospecs.n_pressurewindow = 1.458;
% %         geospecs.n_air = 1.00;
% %         geospecs.n_pressurewindow = 1.52;
% 
%         % jar dimensions in cm
%         geospecs.jar_cylthick = .4; % thickness of cylinder wall
%         geospecs.jar_axthick = .4; % thickness of sphere wall at apex
%         geospecs.jar_cylrad = 15; % outer radius of cylinder
%         geospecs.jar_axrad = 15; % outer radius of sphere (along cylinder axis)
% 
%         geospecs.jar_cyllength = 76;%78.96;
%         geospecs.jar_axrad_top = 15;
%         geospecs.jar_axthick_top = .4;
%         geospecs.jar_bellowsrad = 14.9;
% 
%         % cf3i volume
%         geospecs.target_mass = 38000 + 38000*toplevel;
%         geospecs.target_density = 2;
% 
%         % pressure vessel body (all dimensions inside-dimensions)
%         geospecs.pv_cylbottom = -2.54*13 + sqrt(1-(3.03/(12-(3/16)))^2);%2.54*(2+11+39.11-65.5);%5-15.2*2.54;%14-15;%1*2.54;
%         geospecs.pv_cyllength = 65.5*2.54;%3.5*2.54;
%         geospecs.pv_cylrad = (12-(3/16))*2.54;
%         geospecs.pv_axrad_top = .1;%8.625*2.54;
%         geospecs.pv_axrad_bot = 1;%(2.54*(66.34-1.218)-159.8)/sqrt(1-(3.03/(12-1.218))^2);%8.625*2.54;
%         geospecs.pv_portrad_top = 11.8*2.54;
%         geospecs.pv_portrad_bot = 3.03*2.54;
%         geospecs.pv_top = 2.54*(2+11+39.11); % height relative to jar center, from CM
%         geospecs.pv_bot = 2.54*(2+11+39.11-65.5-4.5); % pipe length to flange surface
% % pv_botreflector = pv_cylbottom - pv_axrad_bot*sqrt(1-(pv_portrad_bot/pv_cylrad)^2);
% 
%         % pressure vessel view ports
%         geospecs.vp_outerrad = .5*6.5*2.54 + tan(pi/6)*10;%.5*9.5*2.54;%5*2.54;
%         geospecs.vp_innerrad = .5*6.5*2.54;%3*2.54;
%         geospecs.vp_winrad = .5*7*2.54;%2.9*2.54;
%         geospecs.vp_conelength = 10;%1.5*tan(pi/3)*2.54;%(2.52-2.5*.5)*2.54;
%         geospecs.vp_innerlength = 2.54*.06;%3*tan(pi/3)*2.54 - 4.4*2.54;%.3*2.54;
%         geospecs.vp_winthick = 2*2.54;
%         geospecs.vp_totallength = 15*2.54;
%         geospecs.vp_height = 2*2.54+11*2.54*toplevel;%+11*2.54;%.75*2.54;
%         geospecs.vp_phi = 60*pi/180;
%         geospecs.vp_lightring_innerrad = 2.54;%1.6*tan(.4/.5);%2.54;
%         geospecs.vp_lightring_outerrad = 2*2.54;
% 
%         % reflector wall
%         geospecs.tworeflectors = false;
%         geospecs.ref_offaxis = 0;%5*2.54;
%         geospecs.ref_cylrad = 25.3;%10*2.54;%8.25*2.54;
%         geospecs.ref_slope_top = 10;%2;%8.25*2.54;
%         geospecs.ref_slope_bot = 2;%8.25*2.54;
%         geospecs.ref_azwidth = 4.23;%4*pi/3;
%         geospecs.ref_cyllength = 2.54*(39.114+11+2+2.75-21.65+0);%3.5*2.54;
%         geospecs.ref_cylbottom = -2.75*2.54;%-1*2.54;
%         geospecs.ref_toplength = .1;%3*2.54;
%         geospecs.ref_botlength = 6.5*2.54;
%         geospecs.ref_slope_bot2 = 1;%8.25*2.54;
%         geospecs.ref_bot2length = 3.75*2.54;%14*2.54;
% 
%         % cameras, position relative to air-side center of viewport, +y towards
%         % chamber
% %updated        
%         geospecs.cam_x = 1.7465;
%         geospecs.cam_y = -41.4691;
%         geospecs.cam_z = 4.7431;
%         geospecs.cam_f = 0.5836;%0.40056;%0.65488; % 23FM65 Tamron
%         geospecs.cam_barreld = [0 0];%[.090877 .51087];%[.1527 .2898];%0.1645; %23FM65 Tamron
%         geospecs.cam_lenstype = 'theta';
%         geospecs.cam_sensorsize = [1088 2048] * .00055;%[491 656] * .00099;%[.44649 .64944];
%         geospecs.cam_resolution = [1088 2048];
% 
%         geospecs.cam_pitch = 0;
% %         geospecs.cam_yaw = 0;
% %         geospecs.cam_roll = 0;


%control / and control t
%% overwrite allparams with fitparams for those params that we are fitting
allparams=zeros(length(which_fit),1);
allparams(which_fit) = fitparams;


%% parse allparams

 geospecs = COUPP60_Geometry_Parms_Jin();

if (which_fit(1)==1)
    geospecs.cam_x = allparams(1);
end
if (which_fit(2)==1)
    geospecs.cam_y = allparams(2);
end
if (which_fit(3)==1)
    geospecs.cam_z = allparams(3);
end
if (which_fit(4)==1)
    geospecs.cam_f = allparams(4);
end
if (which_fit(5)==1)
    geospecs.cam_pitch = allparams(5);
end
if (which_fit(6)==1)
    geospecs.cam_yaw = allparams(6);
end
if (which_fit(7)==1)
    geospecs.cam_roll = allparams(7);
end


%% Create Geometry
[surface_list rays ray_startingpoints pixels] = CreateNew2LGeometry_Jin(geospecs);

%% pixel initialization
pixels_list=zeros(length(rays_list),1);
for ie=1:length(rays_list)
    pixel_ix = find(pixels{1}(:,1)==rays_list(ie,1) & pixels{1}(:,2)==(rays_list(ie,2)+160), 1, 'first');
    if ~isempty(pixel_ix)
        pixels_list(ie)=pixel_ix;
    else
        disp('We have a problem.');
    end
end


%% Create initial rays to trace
inc_rays = rays{1}(pixels_list,:);
ray_sp = ray_startingpoints{1}(pixels_list,:);

%% trace the rays
ray_interfaces = RayTracer(ray_sp, inc_rays, surface_list([1:6 22 23 30 31]), 8, 1e-5, [0.01;2]);

%% find chisq for each traced ray (i.e. pixel)
ki=zeros(length(pixels_list),1);%need changes

x=15*cos(allparams(8));
xl=15*cos(allparams(8)+2*pi/3);
xr=15*cos(allparams(8)-2*pi/3);
y=15*sin(allparams(8));
yl=15*sin(allparams(8)+2*pi/3);
yr=15*sin(allparams(8)-2*pi/3);
z=28-.5*pi*15;

cross_pos=[xl,yl,z+11;xl,yl,z;x,y,z+44;x,y,z+33;x,y,z+22;x,y,z+11;x,y,z;xr,yr,z];%test, only center
%collar=[];
%bot=[0,0,-15];

for ie=1:length(pixels_list),
    r_ie = find(ray_interfaces(target.surface(ie)).ray_index==ie,1,'first');
    if ~isempty(r_ie)
        ki(ie)=sum( (ray_interfaces(target.surface(ie)).intersection_point(r_ie,:)-cross_pos(ie,:)).* ...
            (ray_interfaces(target.surface(ie)).intersection_point(r_ie,:)-cross_pos(ie,:))./ ...
            target.error(ie)*target.error(ie), 2);
    else
        disp(sprintf('Pixel %d missed its intended target!',ie));
        ki(ie) = 100;
    end
end;

%ki(ie+1)=ray_interfaces(target.surfaces())

%% find total chisq

chisq=sum(ki);

