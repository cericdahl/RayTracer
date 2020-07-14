import geospecs
import createGeometry
#import GeometryAndLightsVals
import numpy as np
import time
import finalRays

gs = geospecs.geospecs()

#Fix:
# pi
# matrixes and vectors to numpy arrays


###ASSIGN GEOSPEC ATTRIBUTE VALUES###
#Bubble
gs.bubble_present= True
gs.bubble_radius= 1
gs.bubble_position=np.array([0, 0, -5])

###Indices of Refraction###
gs.n_target = 1.17 #n=1.224 for Ar @ 940nm, 90K
gs.n_jar = 1.4512 #SiO2 @ 940nm
gs.n_hydraulic = 1.22 #1.21 @ 940nm, 146K ;  1.237 @ 7eV, 146K, another said 1.515;
gs.n_pressurewindow = 1.7569 #Al2O3 @ 940nm
gs.n_pressurewall = np.inf
gs.n_air = 1.00 

# jar dimensions in cm
# Outer/first incident Jar (but not container)
gs.ojar_thick = .3 # thickness of cylinder wall 
gs.ojar_cylrad = 12 # outer radius of cylinder 
gs.ojar_axrad = 24.29 # outer radius of sphere (along cylinder axis)
gs.ojar_knucklerad = 4 # don't know how to check
gs.ojar_cyllength = 27.8985
gs.ojar_elevation = 0

# Inner Most(but not pointed) Jar
gs.ijar_thick = 0.5 # thickness of cylinder wall
gs.ijar_cylrad = 10.5 # outer radius of cylinder
gs.ijar_axrad = 21.51 # outer radius of sphere (along cylinder axis)
gs.ijar_knucklerad = 11/3 # don't know how to check
gs.ijar_cyllength = 25.5782 
gs.ijar_elevation = -19.4564

# Camera Position/Orientation
gs.vp_theta = 22.5*np.pi/180 #%20*pi/180;%21.8*pi/180; %correct/22.509
gs.vp_focuselev = -6.531
gs.vp_focuslen = 28.271 # distance from convergence point to CF seal

# Camera Dimensions
gs.vp_win_rad = 1.82372 # mpf is 1.73*.5*2.54; % radius of glass (black on circumference)
gs.vp_air_rad = 1.5875 # radius of air-side can (black on circumference)
#doesn't really matter
gs.vp_can_rad = 2.54 #correct although thats not the only shape used anymore
gs.vp_can_wall = 0.1651
#does again matter
gs.vp_flange_rad = 6.985
gs.vp_win_thick = 0.5080 #measured from location of innersurface
gs.vp_nip_top = -0.6805 #location of innersurface of the window %measured from z=0 for entire cameracan (surface of flange)
gs.vp_can_OAL = 17.3101 #length of the camera can

if 1:
	gs.vp_flange_thick =np.array([7.3025, 1.7526, 1.7272, 1.7272, 1.7272]) #viewport length,bottom half of bottom flange, top half of bottom flange, bottom halp top flange, guess
	gs.vp_nip_rad = 6.6153 # radius of hydraulic-side nipple (black on circumference) %radius of inner viewport cylinder surface
else:
	gs.vp_flange_thick = 2.54*np.array([(2.88-2.382), .69, .625, .625, .625])
	gs.vp_nip_rad = 2.54 # radius of hydraulic-side nipple (black on circumference)

# retroreflector cone and cylindrical surface
gs.rd_rad = 12.5 # reflector-diffuser radius
gs.rd_top = 0 #10
gs.rd_bot = -30
gs.rdcone_top = 8
gs.rdcone_toprad = 8*2.54-.375*2.54

# Weird thing on top
gs.rdtopcone_apex = 16
gs.rdtopcone_rad = 12
gs.rdtopcone_bot = 14 

# Pointed Peice at Bottom
gs.rdbotcone_apex = -15.2
gs.rdbotcone_rad = 10
gs.rdbotcone_bot = -20

# Container Jar (the one the camera sees through) 
gs.pv_top = 9.3230
gs.pv_bot = -83.1799
gs.pv_rad = 8*2.54
gs.pv_thick = .375*2.54 
gs.pv_axrad = (3.07)*2.54 #don't know how to find this

#                                            
# this all seems to modify the light itself 
#                                                
# camera position, orientiation relative to viewport, -z is towards
# chamber, +y is towards jar axis
# (up is cos(vp_theta)\hat{z} + sin(vp_theta)\hat{y} )
gs.cam_x = 0 #not sure if this is worth putting in
gs.cam_y = 0 #not sure if this is worth putting in
gs.cam_z = 0.1395+1.1023 #distance to where rays converge from the inner part of the viewport lens
gs.cam_f = 0.42 #basler 4mm c-mount; 
gs.cam_barreld = np.array([0.015888108817219724,0.04648232478103316]) #2nd,4th order terms % basler 4mm c-mount;  
gs.cam_lenstype = 'theta'
gs.cam_sensorsize = np.array([1024, 1280])*4.8e-4 #I belive these are the physicaldimentions of the ccd
gs.deres = .15 #I beleive this modifies is the resolution of the image
gs.cam_resolution = np.round(gs.deres*(np.array([1024, 1280]))) # this impliments the resolution by modifying the actual number of pixels in the camera by the resolution

gs.cam_pitch = 0
gs.cam_yaw = 0
gs.cam_roll = 0

# parameters for lights 
gs.lights_number=5 #now number per camera
gs.lights_height=-8.5
gs.lights_radius=7.5
gs.lights_nrays=60000 #20000 good for .1
gs.lens_angle=(2/3)*np.pi

# Create light and geometry
[surface_list, rays, startingpoints, pixels] = createGeometry.createGeometry(gs)


#%% Run RayTracer2 For The Camera Rays

t = time.time()

[raytracer_output, throwaway1, throwaway2] = RayTracer2(startingpoints[1], rays[1], surface_list, 18, 1e-5, [1e-5, 100])

#%% Find/Plot Camera Ray Incident Angles
#%more extensive commenting in LED equivilent

n_pixels=np.size(rays[1]) #find N number of rays
n_pixels=n_pixels[1]

pixelid = range(n_pixels_1)

found_endpoints = np.zeros(np.size(pixelid),dtype=bool)

final_rays = finalRays.finalRays()  # Check
final_rays.incoming_ray = np.zeros(n_pixels, 10)
final_rays.intersection_point = np.zeros(n_pixels, 3)
final_rays.surface_normal = np.zeros(n_pixels, 3)
final_rays.ray_index = np.zeros(n_pixels, 1)
final_rays.surface_index = np.zeros(n_pixels, 1)

## What is this?
for n in range(np.size(raytracer_output,0,-1)):
    if np.all(found_endpoints):
        break

    [pixel_present, pixel_ix] = np.isin(pixelid, raytracer_output[n].ray_index)
    these_endpoints = np.logical_and(not found_endpoints, pixel_present)

    found_endpoints = np.logical_or(found_endpoints, these_endpoints)
    final_rays.incoming_ray[pixelid[these_endpoints],:] = raytracer_output[n].incoming_ray[pixel_ix[these_endpoints],:]
    final_rays.intersection_point[pixelid[these_endpoints],:] = raytracer_output[n].intersection_point[pixel_ix[these_endpoints],:]
    final_rays.surface_normal[pixelid[these_endpoints],:] = raytracer_output[n].surface_normal[pixel_ix[these_endpoints],:]
    final_rays.ray_index[pixelid[these_endpoints]] = raytracer_output[n].ray_index[pixel_ix[these_endpoints]]
    final_rays.surface_index[pixelid[these_endpoints]] = raytracer_output[n].surface_index[pixel_ix[these_endpoints]]

final_rays.cos_incident = -np.sum(final_rays.incoming_ray[:,0:2] * final_rays.surface_normal, 1)
final_rays.sin_incident = np.sqrt(1 - final_rays.cos_incident**2)
final_rays.theta_incident = np.arccos(final_rays.cos_incident)

# %plot intersection points colorcoded with incident angles
# %figure(1703)
# %clf;
# %scatter3(final_rays.intersection_point(:,1),final_rays.intersection_point(:,2),final_rays.intersection_point(:,3),10,(180/pi)*final_rays.theta_incident)
# %axis equal
# %c = colorbar;
# %c.Label.String = 'Incident Angle (degrees)';
# %xlabel('x')
# %ylabel('y')
# %zlabel('z')
# %zlim([-30,25])

print(time.time() - t)

# %% Run RayTracer2 For The LED's

t = time.time()

[raytracer_output2, throwaway1, throwaway2] = RayTracer2(startingpoints[2], rays[2], surface_list, 18, 1e-5, [1e-5,100])

# %% Find/Plot LED Incident Angles

n_rays=np.size(rays[2]) # %find N number of rays
n_rays=n_rays[1]

led_id=range(n_rays)
found_endpoints = np.zeroes(np.size(led_id),dtype=bool)  #%prepare a row (maybe should be column?) vector boolean array of length N

final_rays2 = finalRays.finalRays() #Check
final_rays2.incoming_ray = np.zeros(n_rays, 10) #... %does each pixel coorespond to one array?
final_rays2.intersection_point = np.zeros(n_rays, 3)
final_rays2.surface_normal = np.zeros(n_rays, 3)
final_rays2.ray_index = np.zeros(n_rays, 1)
final_rays2.surface_index = np.zeros(n_rays, 1)

for n in range(np.size(raytracer_output2,0,-1)): #%cycles through different scattering steps
    if np.all(found_endpoints):
        break

    # %ray_present is the shape of led_id by virtue of ismember()
    [ray_present, led_ix] = isin(led_id, raytracer_output2[n].ray_index) # %led_ix are the locations of the rays corresponding the pixel ids in raytracer_output2(n).ray_index
    these_endpoints = np.logical_and(not found_endpoints, ray_present)  #%if we have not yet found the point and it was just located then true

    found_endpoints = np.logical_or(found_endpoints, these_endpoints) # %points we already found or points we just found
    final_rays2.incoming_ray[led_id[these_endpoints],:] = raytracer_output2[n].incoming_ray[led_ix[these_endpoints],:]
    final_rays2.intersection_point[led_id[these_endpoints],:] = raytracer_output2[n].intersection_point[led_ix[these_endpoints],:]
    final_rays2.surface_normal[led_id[these_endpoints],:] = raytracer_output2[n].surface_normal[led_ix[these_endpoints],:]
    final_rays2.ray_index[led_id[these_endpoints]] = raytracer_output2[n].ray_index[led_ix[these_endpoints]]
    final_rays2.surface_index[led_id[these_endpoints]] = raytracer_output2[n].surface_index[led_ix[these_endpoints]]

final_rays2.cos_incident = -np.sum(final_rays2.incoming_ray[:,0:2] * final_rays2.surface_normal, 1)
final_rays2.sin_incident = np.sqrt(1 - final_rays2.cos_incident**2)
final_rays2.theta_incident = np.arccos(final_rays2.cos_incident)

# %plot intersection points colorcoded with incident angles
# %figure(1702)
# %clf;
# %scatter3(final_rays2.intersection_point(:,1),final_rays2.intersection_point(:,2),final_rays2.intersection_point(:,3),10,(180/pi)*final_rays2.theta_incident)
# %axis equal
# %c = colorbar;
# %c.Label.String = 'Incident Angle (degrees)';
# %xlabel('x')
# %ylabel('y')
# %zlabel('z')
# %zlim([-30,25])

print(time.time() - t)

# %% Determine Camera Ray Intensities By Relating To LED Rays

t = time.time()

# %calculate and print about how many LED rays there will be for each camera
# %ray to give a sense of how much data the code has to work with
lc_ratio=(gs.lights_number*gs.lights_nrays*3)/(gs.cam_resolution[0]*gs.cam_resolution[1])
print('There are about',lc_ratio, 'LED rays for each camera ray')
# %disp('the are about ' + lc_ratio + ' LED rays for each camera ray')

# %create space in the data for the final pixel rays to store information
# %about what intensity they will receive %and from what led rays
final_rays.pixel_intensity=np.zeros(np.size(pixelid), 1)
# %final_rays.leds=cell(length(pixelid),1)

############# STOP PORTING HERE ##############
#You can do anything you want with the final_rays and final_rays2. To achieve the original model, you do what the rest of the code does:
# compare LED rays to camera rays and pixel match! Then plot.
# For other models, you can use them for whatever you want.
# For Testing purposes, the exact values of the variables should be looked at
# So I didn't think that porting the plotting would be valuable.




############# Nothing after this is ported to Python-- it is all still in MATLAB ############
#
## %create booleans tracking rays on surfaces that matter (for
## %plotting)
#on_retrosurface={false(length(pixelid),1), false(length(led_id),1)};
#
#%make copies of some stucture components as arrays because parfor is picky
#%about structs
#fr_ip=final_rays.intersection_point;
#fr2_ip=final_rays2.intersection_point;
#fr_ir=final_rays.incoming_ray;
#fr2_ir=final_rays2.incoming_ray;
#fr_ti=final_rays.theta_incident;
#fr_ti2=final_rays2.theta_incident;
#%fr_pi=final_rays.pixel_intensity;
#addition_storage=cell(length(led_id),1);
#
#%consider each of the surfaces that we care about
#for n=1:length(surface_list)
#    disc=surface_list(n).description;
#
#    %determine if we should receive any light from the given surface
#    if strcmp(disc,'reflector/diffuser') || strcmp(disc,'reflector/diffuser cone')...
#            || strcmp(disc,'reflector/diffuser botcone')
#        toc
#        disp(disc)
#        %find camera and led rays that end on this surface
#        cam_on_surface = final_rays.surface_index == n | final_rays.surface_index == -n ;
#        led_on_surface = final_rays2.surface_index == n | final_rays2.surface_index == -n;
#
#        %add rays to the boolean of rays on valid surfaces
#        on_retrosurface{1} = on_retrosurface{1} | cam_on_surface;
#        on_retrosurface{2} = on_retrosurface{2} | led_on_surface;
#
#        %find camera and led indices for this surface
#        ci_os=find(cam_on_surface)';
#        li_os=find(led_on_surface)';
#
#        %assign each led ray to a camera ray on the surface
#        for l=li_os
#            %find the distances between the led ray and each camera ray
#            %on the surface
#            closest_cam_info=[0 0]; %first entry camera index second distance
#            for c=ci_os
#                %calculate the distance between rays
#                c_posn = fr_ip(c,:);
#                l_posn = fr2_ip(l,:);
#                distance = norm(c_posn-l_posn); % you would think that norm would normalise the vector, but no; it finds the magnitude
#                %see if this distance is smaller
#                if closest_cam_info(2)==0 || distance < closest_cam_info(2)
#                    closest_cam_info=[c distance];
#                end
#            end
#            %ignore this found better way to do it
#            %add the led index to the list of them for the camera its
#            %nearest to
#            %current_leds=final_rays(closest_cam_info(1)).leds
#            %final_rays(closest_cam_info(1)).leds=[current_leds,l]
#
#            %add to the intensity of the pixel whos ray ended nearest
#            %based on its relation to the current led
#            if closest_cam_info(1) ~= 0 && closest_cam_info(2) < 5                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probably should make this a variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                %find intensity of the led ray
#                int=fr2_ir(l,7);
#
#                %find angle between the led and camera rays
#                cam_rdir=fr_ir(closest_cam_info(1),1:3);
#                led_rdir=fr2_ir(l,1:3);
#                ang_bet=acos(dot(cam_rdir,led_rdir)/(norm(cam_rdir)*norm(led_rdir))); %still can't beleive its called norm
#
#                %find incident angle of the camera ray
#                c_ang_inc=fr_ti(closest_cam_info(1));
#                l_ang_inc=fr_ti2(l);
#
#                %calculate the added intensity to the camera ray and
#                %store the index and added value so it can be
#                %summed outside of this picky parfor loop which
#                %doesn't like indicies that aren't being summed over
#                int_addition=int*exp(-(ang_bet^2)/(2*(10)*(pi/180)) - (c_ang_inc^2)/(2*(45)*(pi/180)) - (closest_cam_info(2)^2)/(2*2^2) - (l_ang_inc^2)/(2*(45)*(pi/180))); %norm ((2*pi*10*45*(2*pi/180)^2)^-1)* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probably should make these variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                addition_storage{l}=[closest_cam_info(1) int_addition];
#            end
#        end
#        %add the new peices of intensity to the pixel intensity
#        for i=1:length(addition_storage)
#            addition_info=addition_storage{i};
#            if ~isempty(addition_info)
#                current_intensity=final_rays.pixel_intensity(addition_info(1));
#                final_rays.pixel_intensity(addition_info(1))=current_intensity + addition_info(2);
#            end
#
#        end
#        %note that pixel_intensity will simpily remain zero for any rays that
#        %end on a surface that doesn't give off light
#
#        %set addition_storage blank again
#        addition_storage=cell(length(led_id),1);
#    end
#end
#
#%plot camera rays on surfaces that give off light
#figure(1704)
#clf;
#scatter3(final_rays.intersection_point(on_retrosurface{1},1),...
#    final_rays.intersection_point(on_retrosurface{1},2),...
#    final_rays.intersection_point(on_retrosurface{1},3),10)
#axis equal
#xlabel('x')
#ylabel('y')
#zlabel('z')
#
#%plot LED rays on surfaces that give off light
#figure(1705)
#clf;
#scatter3(final_rays2.intersection_point(on_retrosurface{2},1),...
#    final_rays2.intersection_point(on_retrosurface{2},2),...
#    final_rays2.intersection_point(on_retrosurface{2},3),10)
#axis equal
#xlabel('x')
#ylabel('y')
#zlabel('z')
#
#%plot camera rays color coded with their pixel intensity
#figure(1706)
#clf;
#scatter3(final_rays.intersection_point(:,1),final_rays.intersection_point(:,2),final_rays.intersection_point(:,3),10,final_rays.pixel_intensity(:))
#axis equal
#c = colorbar;
#c.Label.String = 'Pixel Intensity';
#xlabel('x')
#ylabel('y')
#zlabel('z')
#
#%plot a histogram of the amount of pixels at various intensities
#figure(1707)
#clf;
#histogram(final_rays.pixel_intensity(:))
#
#%plot camera rays less than a standard deviation above the average pixel
#%intensity colorcoded with pixel intensity
#stds_viewed=1;
#pixi_notbig=final_rays.pixel_intensity(:) < mean(final_rays.pixel_intensity(:))...
#    + stds_viewed*std(final_rays.pixel_intensity(:));
#figure(1708)
#clf;
#scatter3(final_rays.intersection_point(pixi_notbig,1),final_rays.intersection_point(pixi_notbig,2),final_rays.intersection_point(pixi_notbig,3),10,final_rays.pixel_intensity(pixi_notbig))
#axis equal
#c = colorbar; %'Ticks',[0,100,300,500,1000,3000],...
#         %'TickLabels',{'0','100','300','500','1000','3000'}
#c.Label.String = stds_viewed + 'Standard Deviation of Pixel Intensity';
#xlabel('x')
#ylabel('y')
#zlabel('z')
#
#toc
#
#%% Make a 2d Histogram of the pixel arrays arranged such as they would be in the camera to create an image
#
#tic
#
#ncolors = 1000;
#% binedges = linspace(0, max(final_rays.theta_incident(cf3i_pixelcut(:))*180/pi), ncolors);
#binedges = linspace(0, 0.5*lc_ratio, ncolors);
#cmap = gray(ncolors); %jet, hot, gray
#[dum bins] = histc(final_rays.pixel_intensity, binedges); %%%%%% Creation of incident angle plot
#bins(bins==0) = length(binedges);
#
#resolution = max(pixels{1}, [], 1);
#composite_image = ones([resolution 3]);
#r_image = zeros(resolution);
#g_image = zeros(resolution);
#b_image = zeros(resolution);
#
#new_b_image = b_image;
#new_g_image = g_image;
#new_r_image = r_image;
#new_r_image(new_r_image==1) = .5;
#
#cf3i_pixelcut=true(resolution);
#new_b_image(cf3i_pixelcut(:)) = cmap(bins(cf3i_pixelcut(:)),3);
#new_g_image(cf3i_pixelcut(:)) = cmap(bins(cf3i_pixelcut(:)),2);
#new_r_image(cf3i_pixelcut(:)) = cmap(bins(cf3i_pixelcut(:)),1);
#
#first_retro_surface = 0;
#miss_retroreflector = abs(final_rays.surface_index)<first_retro_surface;
#
#new_b_image(cf3i_pixelcut(:) & miss_retroreflector(:)) = 0;
#new_g_image(cf3i_pixelcut(:) & miss_retroreflector(:)) = 0;
#new_r_image(cf3i_pixelcut(:) & miss_retroreflector(:)) = 0;
#
#composite_image(:,:,1) = (new_r_image);
#composite_image(:,:,2) = (new_g_image);
#composite_image(:,:,3) = (new_b_image);
#figure(1709);
#clf;
#image(permute(composite_image,[2 1 3]));
#axis equal
#
#colormap(cmap);
#caxis(binedges([1 end]));
#cb = colorbar;
#ch = get(cb,'children');
#set(cb,'ylim',binedges([1 end]));
#set(ch,'Ydata',binedges([1 end]));
#set(cb,'fontsize',16);
#ylabel(cb,'Pixel Intensity');
#
#toc
#








