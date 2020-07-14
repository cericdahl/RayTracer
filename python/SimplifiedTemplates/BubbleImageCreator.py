import geospecs
import createGeometry
# import GeometryAndLightsVals	# What is this?
import numpy as np

gs = geospecs.geospecs()

#Fix:
# pi
# matrixes and vectors to numpy arrays


###ASSIGN GEOSPEC ATTRIBUTE VALUES###
#Bubble
gs.bubble_present = True
gs.bubble_radius = 1
gs.bubble_position = np.array([0, 0, -5])

###Indices of Refraction###
gs.n_target = 1.17 #n=1.224 for Ar @ 940nm, 90K
gs.n_jar = 1.4512 #SiO2 @ 940nm
gs.n_hydraulic = 1.22 #1.21 @ 940nm, 146K ;  1.237 @ 7eV, 146K

# , another said 1.515;
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
	gs.vp_flange_thick = np.array([7.3025, 1.7526, 1.7272, 1.7272, 1.7272]) #viewport length,bottom half of bottom flange, top half of bottom flange, bottom halp top flange, guess
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
gs.pv_axrad = 3.07 * 2.54 #don't know how to find this

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
gs.cam_barreld = np.array([0.015888108817219724, 0.04648232478103316]) #2nd,4th order terms % basler 4mm c-mount;
gs.cam_lenstype = 'theta'
gs.cam_sensorsize = np.array([1024, 1280])*4.8e-4 #I belive these are the physicaldimentions of the ccd
gs.deres = .15 #I believe this modifies the resolution of the image
gs.cam_resolution = round(geospecs.deres*(np.array([1024, 1280]))) # this impliments the resolution by modifying the actual number of pixels in the camera by the resolution

gs.cam_pitch = 0
gs.cam_yaw = 0
gs.cam_roll = 0

# parameters for lights 
gs.lights_number = 5 #now number per camera
gs.lights_height = -8.5
gs.lights_radius = 7.5
gs.lights_nrays = 60000 #20000 good for .1
gs.lens_angle = (2/3)*np.pi

# Create light and geometry
[surface_list, rays, startingpoints, pixels] = createGeometry.createGeometry(gs)

#run RayTracer2
[raytracer_output, throwaway1, throwaway2] = RayTracer2(startingpoints[1], rays[1], surface_list, 18, 1e-5, [1e-5, 100])

####do other calculations and then plot the results here####

print(raytracer_output)
print(throwaway1)
print(throwaway2)








