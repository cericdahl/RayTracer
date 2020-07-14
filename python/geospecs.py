import numpy as np

class geospecs(object):

    def __init__(self): #add deres variable
        #Bubble
        self.bubble_present=False
        self.bubble_radius=1
        self.bubble_position=np.zeros((3), dtype=int)
        
        #Indices of refraction
        self.n_target = 1.17 #% n=1.224 for Ar @ 940nm, 90K
        self.n_jar = 1.4512 #% SiO2 @ 940nm
        self.n_hydraulic = 1.22 #% 1.21 @ 940nm, 146K ;  1.237 @ 7eV, 146K, another said 1.515;
        self.n_pressurewindow = 1.7569 #% Al2O3 @ 940nm
        self.n_pressurewall = np.inf
        self.n_air = 1.00
        
        #Jar dimensions in cm
        #Outer/first incident Jar (but not container)
        self.ojar_thick = .25 #% thickness of cylinder wall
        self.ojar_cylrad = 7.5 #% outer radius of cylinder
        self.ojar_axrad = 15 #% outer radius of sphere (along cylinder axis)
        self.ojar_knucklerad = 2.5
        self.ojar_cyllength = 40
        self.ojar_elevation = 20
        
        #Inner Most(but not pointed) Jar
        self.ijar_thick = .25 #% thickness of cylinder wall
        self.ijar_cylrad = 6.5 #% outer radius of cylinder
        self.ijar_axrad = 13 #% outer radius of sphere (along cylinder axis)
        self.ijar_knucklerad = 2.5
        self.ijar_cyllength = 20
        self.ijar_elevation = 0
        
        #Camera Position/Orientation
        self.vp_theta = 22*np.pi/180 #%21.8*pi/180;
        self.vp_focuselev = -10
        self.vp_focuslen = 33.5 #% distance from convergence point to CF seal
        
        #Camera Dimentions
        self.vp_win_rad = 1.73*.5*2.54 #% radius of glass (black on circumference)
        self.vp_air_rad = 1.25*.5*2.54 #% radius of air-side can (black on circumference)
        #doesn't really matter
        self.vp_can_rad = 2*2.54
        self.vp_can_wall = .125*2.54
        #does again matter
        self.vp_flange_rad = 3.375*2.54
        self.vp_nip_rad = 1.75*.5*2.54 #% radius of hydraulic-side nipple (black on circumference)
        self.vp_win_thick = .25*2.54
        self.vp_nip_top = .5
        self.vp_can_OAL = 6*2.54
        self.vp_flange_thick = np.full((5,), 1)*2.54*.5
        
        #retroreflector cone and cylindrical surface
        self.rd_rad = 12 #% reflector-diffuser radius
        self.rd_top = 100
        self.rd_bot = 0
        self.rdcone_top = 120
        self.rdcone_toprad = 16
        
        #Weird thing on top
        self.rdtopcone_apex = 150
        self.rdtopcone_rad = 10.5
        self.rdtopcone_bot = -20

        #Pointed Peice at Bottom
        self.rdbotcone_apex = -15.2
        self.rdbotcone_rad = 10.5
        self.rdbotcone_bot = -20
        
        #Container Jar (the one the camera sees through)
        self.pv_bot = -20
        self.pv_top = +100
        self.pv_rad = 30
        self.pv_thick = 1
        self.pv_axrad = 15

		##########
		#this all seems to modify the light itself
		##########    
		#camera position, orientiation relative to viewport, -z is towards
		#chamber, +y is towards jar axis
        #(up is cos(vp_theta)\hat{z} + sin(vp_theta)\hat{y} )
        self.cam_x = 0
        self.cam_y = 0
        self.cam_z = 5
        self.cam_f = .8
        self.cam_barreld = 0  #%[-9.21849495e-11  3.32394113e-22  1.15845150e-03 -4.29027882e-17]
        self.cam_lenstype = "theta"
        self.cam_sensorsize = np.full((2,), .1)
        self.cam_resolution = np.array([480, 640])
        self.deres = None ################# NOT IN Create.py
        
        self.cam_pitch = 0
        self.cam_yaw = 0
        self.cam_roll = 0
        
        #parameters for lights
        self.lights_number=5 #%now number per camera
        self.lights_height=-8.5
        self.lights_radius=7.5
        self.lights_nrays=100
        self.lens_angle=(2/3)*np.pi

