# % function [surface_list rays ray_startingpoints pixels] = Create30LGeometry()
# %
# % This function creates a structure array of surfaces to be used by
# % RayTracer.  Follow this architecture to create any geometry you like.
# %
# % Each surface has seven fields.  The first (description) describes the geometry in english.
#   Skipping the first 5 for a second, let's describe what the 6th and 7th
#   parameters are. The 6th parameter, (shape) defines the exact shape out
#   of the following valid options: "plane", "sphere", "cylinder",#
#   "torus". This will be used to call the proper RayToXXXXX when the core
#   engine runs. The 7th parameter (param_list) defines the exact
#   geometric shape of the surface. For example:

# %  a surface with shape of "cylinder" and param_list of [[0, 0, 0], [0, 0, 1], 10]
# % defines a cylinder on the z-axis with radius 10.  See all RayToXXXXX
# % functions for other possible shapes (and
# % create your own if you desire). Right now, for testing purposes, only Cylinder exists
# %
# % The second field (inbounds_function) defines the bounds of the surface,
# % and is a function handle to an anonymous function that inputs an N-by-3-by-M 
# % array and outputs an N-by-M logical.  It can be assumed that all input
# % points are on the surface defined by intersect_function, giving true if
# % the input point is contained within the bounds, false otherwise.  For
# % example, from MATLAB:
# %
# %  @(p)(reshape( ...
# %      (p(:,3,:)>20) & (p(:,3,:)<80) & (atan2(p(:,2,:),p(:,1,:))>0), ...
# %      size(p,1), [] ));
# %
# % would cut the above cylinder in half along the xz plane and truncate it 
# % at 20 and 80 in the z-coordinate.  (Generically, p is an N-by-3-by-M
# % matrix, and the output of the function should be an N-by-M boolean.)
# %
# % The second and third fields are n_outside and n_inside, and give the
# % indices of refraction on the two sides of the surface.  What is 'inside'
# % and 'outside' is defined in the RayToXXXXX function -- for spheres and
# % cylinders, this is obvious, for planes less so.  See the documentation
# % for each RayToXXXXX function for details.  Also, setting n to inf makes
# % that side a perfect conductor (in terms of calculating reflection and
# % polarization, at least).
# %
# % The fourth field is surface type, and may be 'normal', 'diffuse', or
# % 'retro'.  For 'diffuse' and 'retro' surfaces, the normal direction
# % returned by the intersect_function is replaced by a random direction
# % within pi/4 of normal or the reverse of the incoming ray, respectively.
# %
# % The fifth field is an absorption coefficient -- RayTracer will multiply
# % the intensity of both the reflected and refracted rays coming from this
# % surface by 1-absorption.
# %

#import stuff
import numpy as np
import surface
import GenerateRaysFromCamera
import geospecs
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib

#function [surface_list, rays, ray_startingpoints, pixels] = CreateArBCGeometry_WithLights_Split(geospecs)
def createGeometry(gs):
    #I don't know what this does, why would geospecs ever be empty
    #and why do you have to
    # # %% set defaults
    # if nargin<1 || isempty(geospecs) || ~isstruct(geospecs)
    #     geospecs = struct();
    # end


    # #redefined wayy below
    # surface_list = struct( ...
    #     'description', {}, ...
    #     'intersect_function', {}, ...
    #     'inbounds_function', {}, ...
    #     'n_outside', {}, ...
    #     'n_inside', {}, ...
    #     'surface_type', {}, ...
    #     'absorption', {});

    #Instantiate np.ndarray)() for the outputs: rays, ray_startingpoints, pixels
    # rays = cell(1,2);
    # ray_startingpoints = cell(1,2);
    # pixels = cell(1,2);
    #cell is like a list of ndarrays
    rays = [] #np.empty([1,2],dtype=object)
    ray_startingpoints = [] #np.empty([1,2],dtype=object)
    pixels = [] #np.empty([1,2],dtype=object)


    # %% useful equations
    # % n0 = index of refraction at known density
    # % r = rho/rho0 (rho = density now, rho0 = density corresponding to n0)
    # % returns index of refraction at density rho
    #clausius_mossotti = @(n0, r)(sqrt(((1 + 2*r).*n0.*n0 + 2 - 2*r)./((1 - r).*n0.*n0 + 2 + r)));



    #All of these values should be moved to be the default values of the geospecs class

    ################# MOVE TO GEOSPECS DEFAULT VALUES #################
    # # %% indices of refraction and dimensions used below
    # n_target = 1.17 #% n=1.224 for Ar @ 940nm, 90K
    # n_jar = 1.4512 #% SiO2 @ 940nm
    # n_hydraulic = 1.22 #% 1.21 @ 940nm, 146K ;  1.237 @ 7eV, 146K, another said 1.515;
    # n_pressurewindow = 1.7569 #% Al2O3 @ 940nm
    # n_pressurewall = np.inf
    # n_air = 1.00

    # # %Bubble
    # bubble_present=False
    # bubble_radius=1
    # bubble_position=np.zeros((3,), dtype=int)

    # # %%%%                                  %%%%
    # # %%%% Actual Camera Dimentions in here %%%%
    # # %%%%                                  %%%%
    # # % jar dimensions in cm
    # ojar_thick = .25 #% thickness of cylinder wall
    # ojar_cylrad = 7.5 #% outer radius of cylinder
    # ojar_axrad = 15 #% outer radius of sphere (along cylinder axis)
    # ojar_knucklerad = 2.5
    # ojar_cyllength = 40
    # ojar_elevation = 20

    # ijar_thick = .25 #% thickness of cylinder wall
    # ijar_cylrad = 6.5 #% outer radius of cylinder
    # ijar_axrad = 13 #% outer radius of sphere (along cylinder axis)
    # ijar_knucklerad = 2.5
    # ijar_cyllength = 20
    # ijar_elevation = 0

    # vp_theta = 22*pi/180 #%21.8*pi/180;
    # vp_focuselev = -10
    # vp_focuslen = 33.5 #% distance from convergence point to CF seal

    # vp_win_rad = 1.73*.5*2.54 #% radius of glass (black on circumference)
    # vp_air_rad = 1.25*.5*2.54 #% radius of air-side can (black on circumference)
    # vp_can_rad = 2*2.54
    # vp_can_wall = .125*2.54
    # vp_flange_rad = 3.375*2.54
    # vp_nip_rad = 1.75*.5*2.54 #% radius of hydraulic-side nipple (black on circumference)
    # vp_win_thick = .25*2.54
    # vp_nip_top = .5
    # vp_can_OAL = 6*2.54
    # vp_flange_thick = np.full((5,), 1)*2.54*.5

    # rd_rad = 12 #% reflector-diffuser radius
    # rd_top = 100
    # rd_bot = 0
    # rdcone_top = 120
    # rdcone_toprad = 16
    # rdtopcone_apex = 150
    # rdtopcone_rad = 10.5
    # rdtopcone_bot = -20
    # rdbotcone_apex = -15.2
    # rdbotcone_rad = 10.5
    # rdbotcone_bot = -20
            
    # pv_bot = -20
    # pv_top = +100
    # pv_rad = 30
    # pv_thick = 1
    # pv_axrad = 15

    # # %%%%                                       %%%%
    # # %%%% Modifies the Light not the camera can %%%%
    # # %%%%                                       %%%%
    # # % camera position, orientiation relative to viewport, -z is towards
    # # % chamber, +y is towards jar axis
    # # % (up is cos(vp_theta)\hat{z} + sin(vp_theta)\hat{y} )
    # cam_x = 0
    # cam_y = 0
    # cam_z = 5 
    # cam_f = .8
    # cam_barreld = 0  #%[-9.21849495e-11  3.32394113e-22  1.15845150e-03 -4.29027882e-17]
    # cam_lenstype = "theta"
    # cam_sensorsize = np.full((2,), .1)
    # cam_resolution = np.array([480, 640])

    # cam_pitch = 0
    # cam_yaw = 0
    # cam_roll = 0

    # # %%% parameters for lights %%%
    # lights_number=5 #%now number per camera
    # lights_height=-8.5
    # lights_radius=7.5
    # lights_nrays=100
    # lens_angle=(2/3)*np.pi
    ################# END MOVE TO GEOSPECS DEFAULT VALUES #################


    # %% apply geospecs

    # #What does this do?
    # It takes all of the geospecs values, sees if any of the fields overlap with the local variables
    
    
    ############### REMOVE THis ################
    # fn = fieldnames(geospecs)
    # for n=1:length(fn)
    #     if ~isempty(geospecs.(fn{n}))
    #         eval([fn{n} '=geospecs.(fn{n});'])
    #     end
    # end

    #turn all attributes of gs into a dict
    # gsDict = gs.__dict__
    # #get keys of dict
    # fn = gsDict.keys()
    # #loop through attributes
    # for field in fn
    #     if gsDict[field]
    #         #What does this do?
    #         #eval([field '=gsDict[field]'])
    #     end
    # end
    ####################################


    # %% derived dimensions
    vp_s = (gs.vp_focuslen - gs.vp_nip_top) * np.sin(gs.vp_theta) #% radial position of air-side center of viewport
    vp_elev = (gs.vp_focuslen - gs.vp_nip_top) * np.cos(gs.vp_theta) + gs.vp_focuselev #% vertical position of air-side center of viewport

    t_o = np.array([0, gs.ojar_thick])
    t_i = np.array([0, gs.ijar_thick])

    #append arrays rather than make multidimensional
    r1 = np.array([gs.ojar_cylrad - t_o, gs.ijar_cylrad - t_i])
    r2 = np.array([gs.ojar_knucklerad - t_o, gs.ijar_knucklerad - t_i])
    r3 = np.array([gs.ojar_axrad - t_o, gs.ijar_axrad - t_i])
    

    s = r3*(r1-r2)/(r3-r2) #% axis to knuckle-dome transition

    # #CHECK: make sure 1-array is correct
    z = r2*np.sqrt(1 - (s/r3)**2) #%  equator to knuckle-dome transition

    d = r3*z* ((1/r3)-(1/r2)) #% equator to dome sphere center

    cam_pixelpitch = gs.cam_sensorsize / gs.cam_resolution #%think these are the physical dimetnions of a given pixel

    vp_axis = np.array([0, -np.sin(gs.vp_theta), np.cos(gs.vp_theta)])
    vp_center = np.array([0, -vp_s, vp_elev])


    head_out_Q = np.array([[gs.pv_rad**(-2), 0, 0], [0, gs.pv_rad**(-2), 0], [0, 0, gs.pv_axrad**(-2)]])
    head_in_Q = np.array([[(gs.pv_rad-gs.pv_thick)**(-2), 0, 0],[0, (gs.pv_rad-gs.pv_thick)**(-2), 0],[0, 0, (gs.pv_axrad-gs.pv_thick)**(-2)]])
    head_out_P = np.array([0, 0, -2*gs.pv_top* (gs.pv_axrad**(-2))])
    head_in_P = np.array([0, 0, -2*gs.pv_top* ((gs.pv_axrad-gs.pv_thick)**(-2))])
    head_out_R = (gs.pv_top / gs.pv_axrad)**2 - 1
    head_in_R = (gs.pv_top / (gs.pv_axrad - gs.pv_thick))**2 - 1

    rd_cone_b = (gs.rdcone_toprad - gs.rd_rad) / (gs.rdcone_top - gs.rd_top)
    rd_cone_z0 = gs.rd_top - (gs.rd_rad / rd_cone_b)
    rd_cone_Q = np.array([[1, 0, 0],[0, 1, 0],[0, 0, -rd_cone_b**2]])
    rd_cone_P = np.array([0, 0, 2*(rd_cone_b**2)*rd_cone_z0])
    rd_cone_R = -(rd_cone_b * rd_cone_z0)**2

    rd_stcone_b = (gs.rdcone_toprad - gs.rdtopcone_rad) / (gs.rdtopcone_bot - gs.rdcone_top)
    rd_stcone_z0 = gs.rdtopcone_bot + (gs.rdtopcone_rad / rd_stcone_b)
    rd_stcone_Q = np.array([[1, 0, 0],[0, 1, 0],[0, 0, -rd_stcone_b**2]])
    rd_stcone_P = np.array([0, 0, 2*(rd_stcone_b**2)*rd_stcone_z0])
    rd_stcone_R = -(rd_stcone_b * rd_stcone_z0)**2

    rd_topcone_b = gs.rdtopcone_rad / (gs.rdtopcone_apex - gs.rdtopcone_bot)
    rd_topcone_Q = np.array([[1, 0, 0],[0, 1, 0],[0, 0, -rd_topcone_b**2]])
    rd_topcone_P = np.array([0, 0, 2*(rd_topcone_b**2)*gs.rdtopcone_apex])
    rd_topcone_R = -(rd_topcone_b * gs.rdtopcone_apex)**2

    rd_botcone_b = gs.rdbotcone_rad / (gs.rdbotcone_apex - gs.rdbotcone_bot)
    rd_botcone_Q = np.array([[1, 0, 0],[0, 1, 0],[0, 0, -rd_botcone_b**2]])
    rd_botcone_P = np.array([0, 0, 2*(rd_botcone_b**2)*gs.rdbotcone_apex])
    rd_botcone_R = -(rd_botcone_b * gs.rdbotcone_apex)**2

    # %% surface list (where the geometry is actually constructed)
    # %% define surface_list dictionary
    # #Now, surfaceList is a list of classes so we can instantiate multiple instances
    #surface_list ={'description': None, 'intersect_function': None, 'inbounds_function': None, 'n_outside':None, 'n_inside':None, 'surface_type':None, 'absorption':None}
    surface_list = []
    
    # %%% first four silica cylinders
    #surface_list['description'] = 'inside surface of inner quartz cylinder'
    insideInnerCyl = surface.surface()
    insideInnerCyl.description ='inside surface of inner quartz cylinder'
    #surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
         # [0 0 0], [0 0 1], r1(4));
    # Is there actually an r1[5]?
    ## if we add this as a member method we need to pass parameters later
    #insideInnerCyl.intersect_function = RayToCylinder(sp,indir, numpy.array([0 0 0]), numpy.array([0 0 1]), r1(4))
    insideInnerCyl.shape = "cylinder"
    insideInnerCyl.param_list = [[0, 0, 0], [0, 0, 1], r1[3]]
    insideInnerCyl.inbounds_function = lambda p:np.reshape((p[:,2,:]<gs.ijar_elevation) and (p[:,2,:]>=(gs.ijar_elevation-gs.ijar_cyllength)), (np.size(p,0),-1))
    insideInnerCyl.n_outside = gs.n_jar
    insideInnerCyl.n_inside = gs.n_hydraulic
    insideInnerCyl.surface_type = 'normal'
    insideInnerCyl.absorption = 0
    surface_list.append(insideInnerCyl)

    outsideInnerCyl = surface.surface()
    outsideInnerCyl.description = 'outside surface of inner quartz cylinder';
    outsideInnerCyl.shape = "cylinder"
    outsideInnerCyl.param_list = [[0, 0, 0], [0, 0, 1], r1[2]]
    outsideInnerCyl.inbounds_function = lambda p:np.reshape((p[:,2,:]<gs.ijar_elevation) and (p[:,2,:]>=(gs.ijar_elevation-gs.ijar_cyllength)), (np.size(p,0),-1))
    outsideInnerCyl.n_outside = gs.n_target
    outsideInnerCyl.n_inside = gs.n_jar
    outsideInnerCyl.surface_type = 'normal'
    outsideInnerCyl.absorption = 0
    surface_list.append(outsideInnerCyl)

    insideOuterCyl = surface.surface()
    insideOuterCyl.description = 'inside surface of outer quartz cylinder'
    insideOuterCyl.shape = 'cylinder'
    insideOuterCyl.param_list = [[0, 0, 0], [0, 0, 1], r1[1]]
    insideOuterCyl.inbounds_function = lambda p:np.reshape((p[:,2,:]<gs.ojar_elevation) and (p[:,2,:]>=(gs.ojar_elevation-gs.ojar_cyllength)), (np.size(p,0),-1));
    insideOuterCyl.n_outside = gs.n_jar
    insideOuterCyl.n_inside = gs.n_target
    insideOuterCyl.surface_type = 'normal'
    insideOuterCyl.absorption = 0
    surface_list.append(insideOuterCyl)

    outsideOuterCyl = surface.surface()
    outsideOuterCyl.description = 'outside surface of outer quartz cylinder'
    outsideOuterCyl.shape = 'cylinder'
    outsideOuterCyl.param_list = [[0, 0, 0], [0, 0, 1], r1[0]]
    outsideOuterCyl.inbounds_function = lambda p: np.reshape((p[:,2,:]<gs.ojar_elevation) and (p[:,2,:]>=(gs.ojar_elevation-gs.ojar_cyllength)), (np.size(p,0),-1));
    outsideOuterCyl.n_outside = gs.n_hydraulic
    outsideOuterCyl.n_inside = gs.n_jar
    outsideOuterCyl.surface_type = 'normal'
    outsideOuterCyl.absorption = 0
    surface_list.append(outsideOuterCyl)

    # %%% then four silica domes
    insideInnerHemi = surface.surface()
    insideInnerHemi.description = 'inside surface of inner quartz hemisphere'
    # surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    #     [0 0 ijar_elevation+d(4)], r3(4));
    insideInnerHemi.shape = 'sphere'
    insideInnerHemi.param_list = [[0, 0, gs.ijar_elevation+d[3]], r3[3]]
    # surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(4)+ijar_elevation)), size(p,1), [] ));
    insideInnerHemi.inbounds_function = lambda p:np.reshape((p[:,2,:]>(z[3]+gs.ijar_elevation)), (np.size(p,0),-1));
    insideInnerHemi.n_outside = gs.n_jar
    insideInnerHemi.n_inside = gs.n_hydraulic
    insideInnerHemi.surface_type = 'normal'
    insideInnerHemi.absorption = 0
    surface_list.append(insideInnerHemi)

    outsideInnerHemi = surface.surface()
    outsideInnerHemi.description = 'outside surface of inner quartz hemisphere'
    # surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
        # [0 0 ijar_elevation+d(3)], r3(3));
    outsideInnerHemi.shape = 'sphere'
    outsideInnerHemi.param_list = [[0, 0, gs.ijar_elevation+d[2]], r3[2]]
    # surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(3)+ijar_elevation)), size(p,1), [] ));
    outsideInnerHemi.inbounds_function = lambda p:np.reshape((p[:,2,:]>(z[2]+gs.ijar_elevation)), (np.size(p,0),-1))
    outsideInnerHemi.n_outside = gs.n_target
    outsideInnerHemi.n_inside = gs.n_jar
    outsideInnerHemi.surface_type = 'normal'
    outsideInnerHemi.absorption = 0
    surface_list.append(outsideInnerHemi)

    insideOuterHemi = surface.surface()
    insideOuterHemi.description = 'inside surface of outer quartz hemisphere'
    # surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    #     [0 0 ojar_elevation+d(2)], r3(2));
    insideOuterHemi.shape = 'sphere'
    insideOuterHemi.param_list = [[0, 0, gs.ojar_elevation+d[1]], r3[1]]
#    surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(2)+ojar_elevation)), size(p,1), [] ));
    insideOuterHemi.inbounds_function = lambda p:np.reshape((p[:,2,:]>(z[1]+gs.ojar_elevation)), (np.size(p,0),-1))
    insideOuterHemi.n_outside = gs.n_jar
    insideOuterHemi.n_inside = gs.n_target
    insideOuterHemi.surface_type = 'normal'
    insideOuterHemi.absorption = 0
    surface_list.append(insideOuterHemi)

    outsideOuterHemi = surface.surface()
    outsideOuterHemi.description = 'outside surface of outer quartz hemisphere'
#    surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
#       ;
    outsideOuterHemi.shape = 'sphere'
    outsideOuterHemi.param_list = [[0, 0, gs.ojar_elevation+d[0]], r3[0]]
#    surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(z(1)+ojar_elevation)), size(p,1), [] ));
    outsideOuterHemi.inbounds_function = lambda p:np.reshape((p[:,2,:]>(z[0]+gs.ojar_elevation)), (np.size(p,0),-1))
    outsideOuterHemi.n_outside = gs.n_hydraulic
    outsideOuterHemi.n_inside = gs.n_jar
    outsideOuterHemi.surface_type = 'normal'
    outsideOuterHemi.absorption = 0
    surface_list.append(outsideOuterHemi)

    # %%% now four silica knuckles
    insideInnerKnuckle = surface.surface()
    insideInnerKnuckle.description = 'inside surface of inner quartz knuckle'
#    surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
#        [0 0 ijar_elevation], [0 0 1], r1(4)-r2(4), r2(4));
    insideInnerKnuckle.shape = 'torus'
    insideInnerKnuckle.param_list =[[0, 0, gs.ijar_elevation], [0, 0, 1], r1[3]-r2[3], r2[3]]
#    surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>ijar_elevation & p(:,3,:)<=(z(4)+ijar_elevation) & ...
#        (p(:,1,:).^2+p(:,2,:).^2)>((r1(4)-r2(4))^2) ), size(p,1), [] ));
################## p[:,1,:]^2 may cause issues
    insideInnerKnuckle.inbounds_function = lambda p:np.reshape((p[:,2,:]>gs.ijar_elevation and p[:,2,:]<=(z[3]+gs.ijar_elevation) and (p[:,0,:]^2+p[:,1,:]^2)>((r1[3]-r2[3])^2)), (np.size(p,0),-1))
    insideInnerKnuckle.n_outside = gs.n_jar
    insideInnerKnuckle.n_inside = gs.n_hydraulic
    insideInnerKnuckle.surface_type = 'normal'
    insideInnerKnuckle.absorption = 0
    surface_list.append(insideInnerKnuckle)

    outsideInnerKnuckle = surface.surface()
    outsideInnerKnuckle.description = 'outside surface of inner quartz knuckle'
#    surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
#        [0 0 ijar_elevation], [0 0 1], r1(3)-r2(3), r2(3));
    outsideInnerKnuckle.shape = 'torus'
    outsideInnerKnuckle.param_list = [[0, 0, gs.ijar_elevation], [0, 0, 1], r1[2]-r2[2], r2[2]]
#    surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>ijar_elevation & p(:,3,:)<=(z(3)+ijar_elevation) & ...
#        (p(:,1,:).^2+p(:,2,:).^2)>((r1(3)-r2(3))^2) ), size(p,1), [] ));
    outsideInnerKnuckle.inbounds_function = lambda p:np.reshape((p[:,2,:]>gs.ijar_elevation and p[:,2,:]<=(z[2]+gs.ijar_elevation) and (p[:,0,:]^2+p[:,1,:]^2)>((r1[2]-r2[2])^2) ), (np.size(p,0),-1))
    outsideInnerKnuckle.n_outside = gs.n_target
    outsideInnerKnuckle.n_inside = gs.n_jar
    outsideInnerKnuckle.surface_type = 'normal'
    outsideInnerKnuckle.absorption = 0
    surface_list.append(outsideInnerKnuckle)

    insideOuterKnuckle = surface.surface()
    insideOuterKnuckle.description = 'inside surface of outer quartz knuckle'
#    surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
#        [0 0 ojar_elevation], [0 0 1], r1(2)-r2(2), r2(2));
    insideOuterKnuckle.shape = 'torus'
    insideOuterKnuckle.param_list = [[0, 0, gs.ojar_elevation], [0, 0, 1], r1[1]-r2[1], r2[1]]
#    surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>ojar_elevation & p(:,3,:)<=(z(2)+ojar_elevation) & ...
#        (p(:,1,:).^2+p(:,2,:).^2)>((r1(2)-r2(2))^2) ), size(p,1), [] ));
    insideOuterKnuckle.inbounds_function = lambda p:np.reshape((p[:,2,:]> gs.ojar_elevation and p[:,2,:]<=(z[1]+gs.ojar_elevation) and (p[:,0,:]^2+p[:,1,:]^2)>((r1[1]-r2[1])^2)), (np.size(p,0),-1))
    insideOuterKnuckle.n_outside = gs.n_jar
    insideOuterKnuckle.n_inside = gs.n_target
    insideOuterKnuckle.surface_type = 'normal'
    insideOuterKnuckle.absorption = 0
    surface_list.append(insideOuterKnuckle)
    
    outsideOuterKnuckle = surface.surface()
    outsideOuterKnuckle.description = 'outside surface of outer quartz knuckle'
#    surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
#        [0 0 ojar_elevation], [0 0 1], r1(1)-r2(1), r2(1));
    outsideOuterKnuckle.shape = 'torus'
    outsideOuterKnuckle.param_list = [[0, 0, gs.ojar_elevation], [0, 0, 1], r1[0]-r2[0], r2[0]]
    outsideOuterKnuckle.inbounds_function = lambda p:np.reshape((p[:,2,:]>gs.ojar_elevation and p[:,2,:]<=(z[0]+gs.ojar_elevation) and (p[:,0,:]^2+p[:,1,:]^2)>((r1[0]-r2[0])^2) ), (np.size(p,0),-1))
    outsideOuterKnuckle.n_outside = gs.n_hydraulic
    outsideOuterKnuckle.n_inside = gs.n_jar
    outsideOuterKnuckle.surface_type = 'normal'
    outsideOuterKnuckle.absorption = 0
    surface_list.append(outsideOuterKnuckle)

    # %%%% This Mostly Seems to be the Camera %%%%
    # %%% now the viewport
    sightGlassWall = surface.surface()
    sightGlassWall.description = 'sight glass wall'
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        vp_center, vp_axis, vp_air_rad);
    sightGlassWall.shape = 'cylinder'
    sightGlassWall.param_list = [vp_center, vp_axis, gs.vp_air_rad]
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > 0) & ...
#    (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_nip_top + vp_flange_thick(2))), size(p,1), [] ));
    sightGlassWall.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)), 2]) > 0) and (np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)), 2]) <= (gs.vp_nip_top + gs.vp_flange_thick[1])), (np.size(p,0),-1))
    sightGlassWall.n_outside = gs.n_pressurewall
    sightGlassWall.n_inside = gs.n_air
    sightGlassWall.surface_type = 'normal'
    sightGlassWall.absorption = 1
    surface_list.append(sightGlassWall)

    canInnerWall = surface.surface()
    canInnerWall.description = 'camera can inner wall'
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        vp_center, vp_axis, vp_can_rad - vp_can_wall);
    canInnerWall.shape = 'cylinder'
    canInnerWall.param_list = [vp_center, vp_axis, gs.vp_can_rad - gs.vp_can_wall]
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (vp_nip_top + vp_flange_thick(2))) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_can_OAL + vp_flange_thick(2) + vp_nip_top)), size(p,1), [] ));
#    surface_list(end).n_outside = n_pressurewall;
    canInnerWall.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,size(p,2)), 2]) > (gs.vp_nip_top + gs.vp_flange_thick[1])) and (np.sum((p - np.matlib.repmat(vp_center,size(p,0),1,size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,size(p,2)), 2) <= (gs.vp_can_OAL + gs.vp_flange_thick[1] + gs.vp_nip_top)), (np.size(p,0),-1))
    canInnerWall.n_outside = gs.n_pressurewall
    canInnerWall.n_inside = gs.n_air
    canInnerWall.surface_type = 'normal'
    canInnerWall.absorption = 1
    surface_list.append(canInnerWall)
    
    canOuterWall = surface.surface()
    canOuterWall.description = 'camera can outer wall'
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        vp_center, vp_axis, vp_can_rad);
    canOuterWall.shape = 'cylinder'
    canOuterWall.param_list = [vp_center, vp_axis, gs.vp_can_rad]
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3))) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_can_OAL + vp_nip_top + vp_flange_thick(2) - vp_flange_thick(4))), size(p,1), [] ));
    canOuterWall.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)), 2]) > (gs.vp_nip_top + gs.vp_flange_thick[1] + gs.vp_flange_thick[2])) and (np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,size(p,2)), 2]) <= (gs.vp_can_OAL + gs.vp_nip_top + gs.vp_flange_thick[1] - gs.vp_flange_thick[3])), (np.size(p,0),-1))
    canOuterWall.n_outside = 1
    canOuterWall.n_inside = gs.n_pressurewall
    canOuterWall.surface_type = 'normal'
    canOuterWall.absorption = 1
    
    flangeOuterEdge = surface.surface()
    flangeOuterEdge.description = 'flange outer edge'
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        vp_center, vp_axis, vp_flange_rad);
    flangeOuterEdge.shape = 'cylinder'
    flangeOuterEdge.param_list = [vp_center, vp_axis, gs.vp_flange_rad]
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (-vp_flange_thick(1)+vp_nip_top)) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)))) | ...
#        ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4))) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= (vp_nip_top + vp_flange_thick(2) + vp_can_OAL + vp_flange_thick(5)))), size(p,1), [] ));
    flangeOuterEdge.inbounds_function = lambda p:np.reshape(((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2]) > (-gs.vp_flange_thick[0]+gs.vp_nip_top)) and (np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2]) <= (gs.vp_nip_top + gs.vp_flange_thick[1] + gs.vp_flange_thick[2]))) or ((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2]) > (gs.vp_nip_top + gs.vp_flange_thick[1] + gs.vp_can_OAL - gs.vp_flange_thick[3])) and (np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2]) <= (gs.vp_nip_top + gs.vp_flange_thick[1] + gs.vp_can_OAL + gs.vp_flange_thick[4]))),(np.size(p,0),-1))
    flangeOuterEdge.n_outside = 1
    flangeOuterEdge.n_inside = gs.n_pressurewall
    flangeOuterEdge.surface_type = 'normal'
    flangeOuterEdge.absorption = 1
    surface_list.append(flangeOuterEdge)

    windowWall = surface.surface()
    windowWall.description = 'window wall'
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        vp_center, vp_axis, vp_win_rad);
    windowWall.shape = 'cylinder'
    windowWall.param_list = [vp_center, vp_axis, gs.vp_win_rad]
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= 0) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (-vp_win_thick)), ...
#        size(p,1), [] ));
    windowWall.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)), 2]) <= 0) and (np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)), 2]) > (-gs.vp_win_thick)), (np.size(p,0),-1))
    windowWall.n_outside = gs.n_hydraulic
    windowWall.n_inside = gs.n_pressurewindow
    windowWall.surface_type = 'normal'
    windowWall.absorption = 1
    surface_list.append(windowWall)
    
    windowRetainerOuterWall = surface.surface()
    windowRetainerOuterWall.description = 'window retainer outer wall'
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        vp_center, vp_axis, vp_win_rad);
    windowRetainerOuterWall.shape = 'cylinder'
    windowRetainerOuterWall.param_list = [vp_center, vp_axis, gs.vp_win_rad]
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= vp_nip_top) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > 0), ...
#        size(p,1), [] ));
    windowRetainerOuterWall.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2]) <= gs.vp_nip_top) and (np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2]) > 0),(np.size(p,0),-1))
    windowRetainerOuterWall.n_outside = gs.n_hydraulic
    windowRetainerOuterWall.n_inside = gs.n_pressurewall
    windowRetainerOuterWall.surface_type = 'normal'
    windowRetainerOuterWall.absorption = 1
    surface_list.append(windowRetainerOuterWall)

    pressureVesselNippleWall = surface.surface()
    pressureVesselNippleWall.description = 'pressure vessel nipple wall'
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        vp_center, vp_axis, vp_nip_rad);
    pressureVesselNippleWall.shape = 'cylinder'
    pressureVesselNippleWall.param_list=[vp_center, vp_axis, gs.vp_nip_rad]
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) <= vp_nip_top) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)), 2) > (-vp_flange_thick(1)+vp_nip_top)), ...
#        size(p,1), [] ));
    pressureVesselNippleWall.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))*np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2]) <= gs.vp_nip_top) and (np.sum((p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))* np.matlib.repmat(vp_axis,np.size(p,0),1,np.size(p,2)),2) > (-gs.vp_flange_thick[0]+gs.vp_nip_top)), (np.size(p,0),-1))
    pressureVesselNippleWall.n_outside = gs.n_pressurewall
    pressureVesselNippleWall.n_inside = gs.n_hydraulic
    pressureVesselNippleWall.surface_type = 'normal'
    pressureVesselNippleWall.absorption = 1
    surface_list.append(pressureVesselNippleWall)

    #LOOK AT THIS
    viewportAirSide = surface.surface()
    viewportAirSide.description = 'air side of viewport'
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center, vp_axis);
    viewportAirSide.shape = 'plane'
    viewportAirSide.param_list = [vp_center, vp_axis] # 3 element ndarray coordinates
#    surface_list(end).inbounds_function = @(p)(reshape( sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) <= (vp_air_rad^2), size(p,1), [] ));
    viewportAirSide.inbounds_function = lambda p:np.reshape(np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))**2,2]) <= (gs.vp_air_rad**2), (np.size(p,0),-1))  # Anon function inputs a list of pts p and outputs for each a t or f depending on if the point is in part of the shape that's real or not; disks in this case | p[:,1] for x axis
    viewportAirSide.n_outside = gs.n_air
    viewportAirSide.n_inside = gs.n_pressurewindow
    viewportAirSide.surface_type = 'normal'
    viewportAirSide.absorption = 0
    surface_list.append(viewportAirSide)

    viewportHydraulicSide = surface.surface()
    viewportHydraulicSide.description = 'hydraulic side of viewport'
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center - vp_axis*vp_win_thick, vp_axis);
    viewportHydraulicSide.shape = 'plane'
    viewportHydraulicSide.param_list = [vp_center - vp_axis*gs.vp_win_thick, vp_axis]
#    surface_list(end).inbounds_function = @(p)(reshape( sum((p - repmat(vp_center - vp_axis*vp_win_thick,size(p,1),1,size(p,3))).^2,2) <= (vp_win_rad^2), size(p,1), [] ));
    viewportHydraulicSide.inbounds_function = lambda p:np.reshape(np.sum([(p - np.matlib.repmat(vp_center - vp_axis*gs.vp_win_thick,np.size(p,0),1,np.size(p,2)))**2,2]) <= (gs.vp_win_rad**2), (np.size(p,0),-1))
    viewportHydraulicSide.n_outside = gs.n_pressurewindow
    viewportHydraulicSide.n_inside = gs.n_hydraulic
    viewportHydraulicSide.surface_type = 'normal'
    viewportHydraulicSide.absorption = 0
    surface_list.append(viewportHydraulicSide)

    viewportRetainer = surface.surface()
    viewportRetainer.description = 'viewport retainer'
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center, vp_axis);
    viewportRetainer.shape = 'plane'
    viewportRetainer.param_list = [vp_center, vp_axis]
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) > (vp_air_rad^2)) & ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) <= (vp_win_rad^2)), size(p,1), [] ));
    viewportRetainer.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))**2,2]) > (gs.vp_air_rad**2)) and (np.sum([(p - np.matlib.repmat(vp_center,np.size(p,0),1,np.size(p,2)))**2,2]) <= (gs.vp_win_rad**2)), (np.size(p,0),-1))
    viewportRetainer.n_outside = gs.n_pressurewall
    viewportRetainer.n_inside = gs.n_pressurewindow
    viewportRetainer.surface_type = 'normal'
    viewportRetainer.absorption = 1
    surface_list.append(viewportRetainer)

    nippleBottom = surface.surface()
    nippleBottom.description = 'nipple bottom'
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center - vp_axis*(vp_flange_thick(1)-vp_nip_top), vp_axis);
    nippleBottom.shape = 'plane'
    nippleBottom.param_list = [vp_center - vp_axis*(gs.vp_flange_thick[0]-gs.vp_nip_top), vp_axis]
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        (sum((p - repmat(vp_center - vp_axis*(vp_flange_thick(1)-vp_nip_top),size(p,1),1,size(p,3))).^2,2) > (vp_nip_rad^2)) & ...
#        (sum((p - repmat(vp_center - vp_axis*(vp_flange_thick(1)-vp_nip_top),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
    nippleBottom.inbounds_function = lambda p:np.reshape((np.sum([(p - np.matlib.repmat(vp_center - vp_axis*(gs.vp_flange_thick[0]-gs.vp_nip_top),np.size(p,0),1,np.size(p,2)))**2,2]) > (gs.vp_nip_rad**2)) and (np.sum([(p - np.matlib.repmat(vp_center - vp_axis*(gs.vp_flange_thick[0]-gs.vp_nip_top),np.size(p,0),1,np.size(p,2)))**2,2]) <= (gs.vp_flange_rad**2)),(np.size(p,0),-1))
    nippleBottom.n_outside = gs.n_pressurewall
    nippleBottom.n_inside = gs.n_hydraulic
    nippleBottom.surface_type = 'normal'
    nippleBottom.absorption = 1
    surface_list.append(nippleBottom)
    
    
    ######## STOP CONVERTING SURFACES HERE ########


#    surface_list(end+1).description = 'nipple top';
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center + vp_axis*vp_nip_top, vp_axis);
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        (sum((p - repmat(vp_center + vp_axis*vp_nip_top,size(p,1),1,size(p,3))).^2,2) > (vp_win_rad^2)) & ...
#        (sum((p - repmat(vp_center + vp_axis*vp_nip_top,size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
#    surface_list(end).n_outside = n_pressurewall;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'can bot';
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2)), vp_axis);
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2)),size(p,1),1,size(p,3))).^2,2) > (vp_air_rad^2)) & ...
#        (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
#    surface_list(end).n_outside = n_air;
#    surface_list(end).n_inside = n_pressurewall;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'can bot_top';
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)), vp_axis);
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)),size(p,1),1,size(p,3))).^2,2) > (vp_can_rad^2)) & ...
#        (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_flange_thick(3)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
#    surface_list(end).n_outside = 1;
#    surface_list(end).n_inside = n_pressurewall;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'can top_bot';
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4)), vp_axis);
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4)),size(p,1),1,size(p,3))).^2,2) > (vp_can_rad^2)) & ...
#        (sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL - vp_flange_thick(4)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2)), size(p,1), [] ));
#    surface_list(end).n_outside = n_pressurewall;
#    surface_list(end).n_inside = 1;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'can top';
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL), vp_axis);
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2), size(p,1), [] ));
#    surface_list(end).n_outside = n_pressurewall;
#    surface_list(end).n_inside = n_air;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'can very top';
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL + vp_flange_thick(5)), vp_axis);
#    surface_list(end).inbounds_function = @(p)( reshape(...
#        sum((p - repmat(vp_center + vp_axis*(vp_nip_top + vp_flange_thick(2) + vp_can_OAL + vp_flange_thick(5)),size(p,1),1,size(p,3))).^2,2) <= (vp_flange_rad^2), size(p,1), [] ));
#    surface_list(end).n_outside = 1;
#    surface_list(end).n_inside = n_pressurewall;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %%% now other black surfaces to trap rays
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)>=0) & (p(:,2,:)>=0) & (p(:,3,:)>=(rd_bot+(rd_top-rd_bot)/2)) ),  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#
#    # %%%%%%%%%%%%%%%%%%% REPEAT TO SPLIT UP SURFACES  note original also modified %%%%%%%%%%%%%%%%%%%%%%%%%%%
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)<=0) & (p(:,2,:)>=0) & (p(:,3,:)>=(rd_bot+(rd_top-rd_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)>=0) & (p(:,2,:)<=0) & (p(:,3,:)>=(rd_bot+(rd_top-rd_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)<=0) & (p(:,2,:)<=0) & (p(:,3,:)>=(rd_bot+(rd_top-rd_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %switch from top to bottom
#
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)>=0) & (p(:,2,:)>=0) & (p(:,3,:)<=(rd_bot+(rd_top-rd_bot)/2)) ),  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)<=0) & (p(:,2,:)>=0) & (p(:,3,:)<=(rd_bot+(rd_top-rd_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)>=0) & (p(:,2,:)<=0) & (p(:,3,:)<=(rd_bot+(rd_top-rd_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], rd_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_bot) & ...
#        (p(:,3,:)<=rd_top) & (p(:,1,:)<=0) & (p(:,2,:)<=0) & (p(:,3,:)<=(rd_bot+(rd_top-rd_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %%%%%%%%%%%%%%%%%% REPEAT END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)<=0) & (p(:,1,:)<=0) & (p(:,3,:)<=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %%%%%% REPEAT TO TRY TO SPLIT UP SURFACES AND MAKE CODE RUN FASTER note
#    # %%%%%% that original is also modified
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)>=0) & (p(:,1,:)<=0) & (p(:,3,:)<=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)<=0) & (p(:,1,:)>=0) & (p(:,3,:)<=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)>=0) & (p(:,1,:)>=0) & (p(:,3,:)<=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %%switch from bottom to top
#
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)<=0) & (p(:,1,:)<=0) & (p(:,3,:)>=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)>=0) & (p(:,1,:)<=0) & (p(:,3,:)>=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)<=0) & (p(:,1,:)>=0) & (p(:,3,:)>=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_cone_Q, rd_cone_P, rd_cone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rd_top) & ...
#        (p(:,3,:)<rdcone_top) & (p(:,2,:)>=0) & (p(:,1,:)>=0) & (p(:,3,:)>=(rd_top+(rdcone_top-rd_top)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %%%%%%%% REPEATS END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#    surface_list(end+1).description = 'reflector/diffuser strip cone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_stcone_Q, rd_stcone_P, rd_stcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdcone_top) & ...
#        (p(:,3,:)<rdtopcone_bot) & ...
#        ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) - ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)),2).^2)) > (vp_nip_rad^2))) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser topcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_topcone_Q, rd_topcone_P, rd_topcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdtopcone_bot) & ...
#        (p(:,3,:)<rdtopcone_apex) & ...
#        ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) - ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)),2).^2)) > (vp_nip_rad^2))) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)<=0) & (p(:,2,:)<=0) & (p(:,3,:)<=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ),  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %%%%%%%%%%%%%%%%%%%%%%% REPEAT TO SPLIT UP SURFACES note original also
#    # %%%%%%%%%%%%%%%%%%%%%%% modified
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)>=0) & (p(:,2,:)<=0) & (p(:,3,:)<=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)<=0) & (p(:,2,:)>=0) & (p(:,3,:)<=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ),  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)>=0) & (p(:,2,:)>=0) & (p(:,3,:)<=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # % switch from bottom to top
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)<=0) & (p(:,2,:)<=0) & (p(:,3,:)>=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ),  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)>=0) & (p(:,2,:)<=0) & (p(:,3,:)>=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)<=0) & (p(:,2,:)>=0) & (p(:,3,:)>=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ),  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'reflector/diffuser botcone';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        rd_botcone_Q, rd_botcone_P, rd_botcone_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>rdbotcone_bot) & ...
#        (p(:,3,:)<rdbotcone_apex) & (p(:,1,:)>=0) & (p(:,2,:)>=0) & (p(:,3,:)>=(rdbotcone_bot+(rdbotcone_apex-rdbotcone_bot)/2)) ) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_hydraulic;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    # %%%%%%%%%%%%%%%%%%%%%%%%%% REPEAT END %%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#    # %%%% I think the container (the "jar" with the hole for the camera") %%%%
#    surface_list(end+1).description = 'PV - cylinder outer wall';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], pv_rad);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>pv_bot) & ...
#        (p(:,3,:)<pv_top)) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = 1;
#    surface_list(end).n_inside = n_pressurewall;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'PV - cylinder inner wall';
#    surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
#        [0 0 0], [0 0 1], pv_rad - pv_thick);
#    surface_list(end).inbounds_function = @(p)(reshape( ( ...
#        (p(:,3,:)>pv_bot) & ...
#        (p(:,3,:)<pv_top)) ,  ...
#        size(p,1), []));
#    surface_list(end).n_outside = n_pressurewall;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'PV - outer top';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        head_out_Q, head_out_P, head_out_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        (p(:,3,:) > pv_top) & ...
#        ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) - ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)),2).^2)) > (vp_flange_rad^2)), size(p,1), [] ));
#    surface_list(end).n_outside = 1;
#    surface_list(end).n_inside = n_pressurewall;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'PV - inner top';
#    surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
#        head_in_Q, head_in_P, head_in_R);
#    surface_list(end).inbounds_function = @(p)(reshape( ...
#        (p(:,3,:) > pv_top) & ...
#        ((sum((p - repmat(vp_center,size(p,1),1,size(p,3))).^2,2) - ...
#        (sum((p - repmat(vp_center,size(p,1),1,size(p,3))).*repmat(vp_axis,size(p,1),1,size(p,3)),2).^2)) > (vp_flange_rad^2)), size(p,1), [] ));
#    surface_list(end).n_outside = n_pressurewall;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;
#
#    surface_list(end+1).description = 'PV - bot';
#    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
#        [0, 0, pv_bot], [0, 0, -1]);
#    surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2 + p(:,2,:).^2) <= (pv_rad^2), size(p,1), [] ));
#    surface_list(end).n_outside = n_pressurewall;
#    surface_list(end).n_inside = n_hydraulic;
#    surface_list(end).surface_type = 'normal';
#    surface_list(end).absorption = 1;


###################
# Skip all these surfaces
###################

    if gs.bubble_present:
        bubble = surface.surface()
        bubble.description = 'bubble'
        bubble.shape = 'sphere'
        bubble.param_list = [gs.bubble_position, gs.bubble_radius]
#        surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
#            bubble_position, bubble_radius);
#        surface_list(end).inbounds_function = @(p)(reshape(p(:,3,:)>-500, size(p,1), [] ));
        bubble.inbounds_function = lambda p:np.reshape(p[:,2,:]>-500, (np.size(p,0),-1))                ## What is this?
        bubble.n_outside = gs.n_target
        bubble.n_inside = gs.n_air
        bubble.surface_type = 'normal'
        bubble.absorption = 0
        surface_list.append(bubble)

    # %% create ray lists for Camera
    [raydirections, pixelmap] = GenerateRaysFromCamera.GenerateRaysFromCamera(gs.cam_resolution, cam_pixelpitch, .5*(1+gs.cam_resolution), gs.cam_f, gs.cam_pitch + gs.vp_theta - (np.pi/2), gs.cam_yaw, gs.cam_roll, gs.cam_barreld, gs.cam_lenstype)

#    rays{1} = [raydirections repmat([0 0 1 1 0 0 0],size(raydirections,1),1)];
#    pixels{1} = pixelmap;
    
    rays.append([raydirections, np.matlib.repmat([0, 0, 1, 1, 0, 0, 0],np.size(raydirections,0),1)])
    pixels.append(pixelmap)
    

    ray_startingpoints.append(np.matlib.repmat(vp_center + [gs.cam_x, 0, 0] + gs.cam_z*vp_axis + gs.cam_y*np.cross(vp_axis, [1,0,0]),(np.size(raydirections,0),1))) #don't understand how this acheives what I would think it's acheiving


    # %% create ray lists for lights
    # %determine starting points
    ray_startingpoints.append(np.zeros((gs.lights_nrays*gs.lights_number*3,3))) #*3 because there are three cameras that the lights are arranged around
    # %loop through each camera
#    for c=1:3
    for c in range(1,4):
        c_angle=(2*np.pi/3)*c
        
        # %find the rotation matrix to rotate to angle of the camera
#        rot_mat=[cos(c_angle) -sin(c_angle) 0;...
#            sin(c_angle) cos(c_angle) 0;...
#            0 0 1];
        rot_mat=np.array([[np.cos(c_angle), -np.sin(c_angle), 0], [np.sin(c_angle), np.cos(c_angle), 0],[0, 0, 1]])
        
        # %rotate the camera location and axis to the one desired
        c_point=np.transpose((rot_mat*np.transpose(vp_center))) #point for the view port in question around which the set of leds will be arranged (transpose due to vp_center/point boing row vectors)
        c_axis=np.transpose((rot_mat*np.transpose(vp_axis))) #again transposes (adjoints really but everything's real) because all the variables involved are row vectors, but I wrote the rotation matrix for column vectors
        
        # %loop through each light around the camera
#        for n=1:lights_number
        for n in range(1,gs.lights_number+1):
            l_angle=(2*np.pi/gs.lights_number)*n
            
            y=[0, 1, 0]
            # %define an orthonormal basis in which the camera axis is z
            z_prime=c_axis/np.linalg.norm(c_axis) #normalise c_axis (who decided to make norm the function that finds magnitude?)
            y_prime=y-np.dot(y,z_prime)*z_prime #project y onto z plane
            y_prime=y_prime/np/linalg.norm(y_prime) #normalizing
            x_prime=np.cross(y_prime,z_prime) #take a cross product to get x_prime
            x_prime=x_prime/np.linalg.norm(x_prime) #normalizing (not really neccecary)
            # %convert a rotation around z in the prime basis to the standard basis
            p=[np.transpose(x_prime), np.transpose(y_prime), np.transpose(z_prime)] #matrix to convert form prime basis to standard basis
#            rot_mat3d=p*[cos(l_angle) -sin(l_angle) 0;...
#            sin(l_angle) cos(l_angle) 0;...
#            0 0 1]*inv(p); %rotation around z_prime/c_axis in the standard basis
            
            rot_mat3d=p* np.array([np.cos(l_angle), -np.sin(l_angle), 0], [np.sin(l_angle), np.cos(l_angle), 0], [0, 0, 1])*inv(p) #rotation around z_prime/c_axis in the standard basis
            
            
            # %calculate & normalise the projection of c_point onto the c_axis plane
            rad_vec=c_point-np.dot(c_point,c_axis)*c_axis/(np.linalg.norm(c_axis)**2)
            rad_vec=rad_vec/np.linalg.norm(rad_vec) #normalizing
            # %rotate the projection
            rad_vec=np.transpose(rot_mat3d*np.transpose(rad_vec)) #adjoints for reasons above
            
            #specfy the location of the light starting from the camera point
            light_loc=c_point + gs.lights_height*c_axis/np.linalg.norm(c_axis)+gs.lights_radius*rad_vec
            
            # %add the location of that lights rays to the starting point array
            c_startindex=gs.lights_nrays*gs.lights_number*(c-1)
#            ray_startingpoints{2}((c_startindex+(n-1)*lights_nrays+1):(c_startindex + n*lights_nrays),:)=repmat(light_loc,lights_nrays,1);
            ray_startingpoints[1][(c_startindex+(n-1)*gs.lights_nrays):(c_startindex + n*gs.lights_nrays),:] = np.matlib.repmat(gs.light_loc,gs.lights_nrays,1)

    # %specify ray properties starting with a blank array
    rays.append(np.zeros(gs.lights_nrays*gs.lights_number,10))

    # %create directions within the angular region determined by the lens centered in the direction of the camera
#    for c=1:3
    for c in range(1,4):
        # %more or less repeat procees for starting points for each camera
        # %(probably should have done it all at once, but I thought I'd get
        # %myself confused).
        c_angle=(2*np.pi/3)*c
        
        # %find the rotation matrix to rotate to angle of the camera
        rot_mat=np.array([np.cos(c_angle), -np.sin(c_angle), 0], [np.sin(c_angle), np.cos(c_angle), 0], [0, 0, 1])
        
        # %rotate the camera axis to the one desired and create an orthonormal
        # %basis where -c_axis is in the z direction
        c_axis=np.transpose(rot_mat*np.transpose(vp_axis)) #adjoints see above
        y=[0, 1, 0]
        z_prime=-c_axis/np.linalg.norm(c_axis) #here -c_axis because we want z facing the chamber
        y_prime=y-np.dot(y,z_prime)*z_prime
        y_prime=y_prime/np.linalg.norm(y_prime)
        x_prime=np.cross(y_prime,z_prime)
        x_prime=x_prime/np.linalg.norm(x_prime)

        # %generate ray directions
        thetas=2*np.pi*np.random.randn(lights_nrays*lights_number,1)
        oneminuscosphis=(1-np.cos(gs.lens_angle/2))*np.random.randn(gs.lights_nrays*gs.lights_number,1) #I am randomising 1-cosphi as opposed to phi as its probability remains uniform and varies 0,2 as phi varies 0,pi
        phis = np.arccos(-oneminuscosphis+1) #find phis
        c_raydirs=np.multiply(np.sin(phis),np.cos(thetas))*x_prime+ np.multiply(np.sin(phis),np.sin(thetas))*y_prime + (-oneminuscosphis+1)*z_prime #using some elementwise multiplication
        
        # %place ray directions in the ray array
        c_startindex=gs.lights_nrays*gs.lights_number*(c-1)
        c_endindex=gs.lights_nrays*gs.lights_number*c
#        rays{2}((c_startindex+1):c_endindex,1:3)=c_raydirs;
        rays[1][(c_startindex):c_endindex,0:3] = c_raydirs


    # %randomize polarization properties
#    rays{2}(:,4:6)=rand(lights_nrays*lights_number*3,3) - rand(lights_nrays*lights_number*3,3);
#    rays{2}(:,7:10)=repmat([1 0 0 0],lights_nrays*lights_number*3,1);
#    return %remove this to see a plot of where the lights are
    rays[1][:,3:6]=np.random.randn(gs.lights_nrays*gs.lights_number*3,3) - np.random.randn(gs.lights_nrays*gs.lights_number*3,3)
    rays[1][:,6:10]=np.matlib.repmat([1, 0, 0, 0],gs.lights_nrays*gs.lights_number*3,1)
#   return %remove this to see a plot of where the lights are


    # %% Test Scatter Plot of Light Rays
#    figure(1701);
    fig = matplotlib.pyplot.figure(1701)
#    clf;
    matplotlib.pyplot.clf
#    point2=ray_startingpoints{2}(:,1:3)+rays{2}(:,1:3)
    point2=ray_startingpoints[1][:,0:3]+rays[1][:,0:3]
#    plot3(ray_startingpoints{2}(:,1), ray_startingpoints{2}(:,2), ray_startingpoints{2}(:,3), '*r', point2(:,1), point2(:,2), point2(:,3), '*b')
    ax = Axes3D(fig)
    ax.scatter(ray_startingpoints[1][:,0], ray_startingpoints[1][:,1], ray_startingpoints[1][:,2],marker='o',color='red')
    ax.scatter(point2[:,0],point2[:,1], point2[:,2],marker='o',color='blue')
#    axis equal
    ax.axis('equal')
    matplotlib.pyplot.xlabel('x')
    matplotlib.pyplot.ylabel('y')
    matplotlib.pyplot.zlabel('z')
    
    #Return the value
    return [surface_list, rays, ray_startingpoints, pixels]
