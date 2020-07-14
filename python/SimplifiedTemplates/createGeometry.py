# % function [surface_list rays ray_startingpoints pixels] = Create30LGeometry()
# %
# % This function creates a structure array of surfaces to be used by
# % RayTracer.  Follow this architecture to create any geometry you like.
# %
# % Each surface has six fields.  The first (intersect_function) defines the
# % geometry, and is a function handle to an anonymous function, calling a
# % RayToXXXXX function with the appropriate geometry inputs.  For example, 
# %
# %  @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], 10) 
# %
# % defines a cylinder on the z-axis with radius 10.  See all RayToXXXXX
# % functions in the RayTracing directory for other possible shapes (and
# % create your own if you desire).
# %
# % The second field (inbounds_function) defines the bounds of the surface,
# % and is a function handle to an anonymous function that inputs an N-by-3-by-M 
# % array and outputs an N-by-M logical.  It can be assumed that all input
# % points are on the surface defined by intersect_function, giving true if
# % the input point is contained within the bounds, false otherwise.  For
# % example,
# %
# %  @(p)(reshape( ...
# %      (p(:,3,:)>20) & (p(:,3,:)<80) & (atan2(p(:,2,:),p(:,1,:))>0), ...
# %      size(p,1), [] ));
# %
# % would cut the above cylinder in half along the xz plane and truncate it 
# % at 20 and 80 in the z-coordinate.  (Generically, p is an N-by-3-by-M
# % matrix, and the output of the function should be an N-by-M boolean.)
# %
# % The third and fourth fields are n_outside and n_inside, and give the
# % indices of refraction on the two sides of the surface.  What is 'inside'
# % and 'outside' is defined in the RayToXXXXX function -- for spheres and
# % cylinders, this is obvious, for planes less so.  See the documentation
# % for each RayToXXXXX function for details.  Also, setting n to inf makes
# % that side a perfect conductor (in terms of calculating reflection and
# % polarization, at least).
# %
# % The fifth field is surface type, and may be 'normal', 'diffuse', or
# % 'retro'.  For 'diffuse' and 'retro' surfaces, the normal direction
# % returned by the intersect_function is replaced by a random direction
# % within pi/4 of normal or the reverse of the incoming ray, respectively.
# %
# % The sixth field is an absorption coefficient -- RayTracer will multiply
# % the intensity of both the reflected and refracted rays coming from this
# % surface by 1-absorption.
# %
# % 12/16/09, CED

# import libraries
import numpy as np
import surface
import geospecs
import CreateSomeLightRays  # This doesn't actually exist it's just a placeholder for whatever function you use


# import matplotlib
# from mpl_toolkits.mplot3d import Axes3D

# function [surface_list, rays, ray_startingpoints, pixels] = CreateArBCGeometry_WithLights_Split(geospecs)
def createGeometry(gs):
    # Initialize some containers
    rays = []
    ray_startingpoints = []
    pixels = []  # This will return empty because none of the code in this example defines the pixels

    # some useful equations
    # n0 = index of refraction at known density
    # r = rho/rho0 (rho = density now, rho0 = density corresponding to n0)
    # returns index of refraction at density rho
    # clausius_mossotti = @(n0, r)(sqrt(((1 + 2*r).*n0.*n0 + 2 - 2*r)./((1 - r).*n0.*n0 + 2 + r)));

    # Calculate some dimensions based on the values provided by geospecs (gs). These may be useful for defining some geometry
    vp_s = (gs.vp_focuslen - gs.vp_nip_top) * np.sin(gs.vp_theta)  # % radial position of air-side center of viewport
    vp_elev = (gs.vp_focuslen - gs.vp_nip_top) * np.cos(
        gs.vp_theta) + gs.vp_focuselev  # % vertical position of air-side center of viewport

    t_o = np.array([0, gs.ojar_thick])
    t_i = np.array([0, gs.ijar_thick])

    r1 = np.array([gs.ojar_cylrad - t_o, gs.ijar_cylrad - t_i])
    r2 = np.array([gs.ojar_knucklerad - t_o, gs.ijar_knucklerad - t_i])
    r3 = np.array([gs.ojar_axrad - t_o, gs.ijar_axrad - t_i])

    s = r3 * (r1 - r2) / (r3 - r2)  # % axis to knuckle-dome transition

    z = r2 * np.sqrt(1 - (s / r3) ** 2)  # %  equator to knuckle-dome transition

    d = r3 * z * ((1 / r3) - (1 / r2))  # % equator to dome sphere center

    vp_axis = np.array([0, -np.sin(gs.vp_theta), np.cos(gs.vp_theta)])
    vp_center = np.array([0, -vp_s, vp_elev])

    head_out_Q = np.array([[gs.pv_rad ** (-2), 0, 0], [0, gs.pv_rad ** (-2), 0], [0, 0, gs.pv_axrad ** (-2)]])
    head_in_Q = np.array([[(gs.pv_rad - gs.pv_thick) ** (-2), 0, 0], [0, (gs.pv_rad - gs.pv_thick) ** (-2), 0],
                          [0, 0, (gs.pv_axrad - gs.pv_thick) ** (-2)]])
    head_out_P = np.array([0, 0, -2 * gs.pv_top * (gs.pv_axrad ** (-2))])
    head_in_P = np.array([0, 0, -2 * gs.pv_top * ((gs.pv_axrad - gs.pv_thick) ** (-2))])
    head_out_R = (gs.pv_top / gs.pv_axrad) ** 2 - 1
    head_in_R = (gs.pv_top / (gs.pv_axrad - gs.pv_thick)) ** 2 - 1

    rd_cone_b = (gs.rdcone_toprad - gs.rd_rad) / (gs.rdcone_top - gs.rd_top)
    rd_cone_z0 = gs.rd_top - (gs.rd_rad / rd_cone_b)
    rd_cone_Q = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -rd_cone_b ** 2]])
    rd_cone_P = np.array([0, 0, 2 * (rd_cone_b ** 2) * rd_cone_z0])
    rd_cone_R = -(rd_cone_b * rd_cone_z0) ** 2

    rd_stcone_b = (gs.rdcone_toprad - gs.rdtopcone_rad) / (gs.rdtopcone_bot - gs.rdcone_top)
    rd_stcone_z0 = gs.rdtopcone_bot + (gs.rdtopcone_rad / rd_stcone_b)
    rd_stcone_Q = np.array([1, 0, 0], [0, 1, 0], [0, 0, -rd_stcone_b ** 2])
    rd_stcone_P = np.array([0, 0, 2 * (rd_stcone_b ** 2) * rd_stcone_z0])
    rd_stcone_R = -(rd_stcone_b * rd_stcone_z0) ** 2

    rd_topcone_b = gs.rdtopcone_rad / (gs.rdtopcone_apex - gs.rdtopcone_bot)
    rd_topcone_Q = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -rd_topcone_b ** 2]])
    rd_topcone_P = np.array([0, 0, 2 * (rd_topcone_b ** 2) * gs.rdtopcone_apex])
    rd_topcone_R = -(rd_topcone_b * gs.rdtopcone_apex) ** 2

    rd_botcone_b = gs.rdbotcone_rad / (gs.rdbotcone_apex - gs.rdbotcone_bot)
    rd_botcone_Q = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -rd_botcone_b ** 2]])
    rd_botcone_P = np.array([0, 0, 2 * (rd_botcone_b ** 2) * gs.rdbotcone_apex])
    rd_botcone_R = -(rd_botcone_b * gs.rdbotcone_apex) ** 2

    # Define a surface list (where the geometry is actually constructed)
    # surface_list ={'description': None, 'intersect_function': None, 'inbounds_function': None, 'n_outside':None, 'n_inside':None, 'surface_type':None, 'absorption':None}
    surface_list = []

    # Define some geometry (Do as many of these as necessary to define all geometry. Right now, we'll just do 2, + 1 bubble)

    # Define Surface #1
    # instantiate the surface
    insideInnerCyl = surface.surface()
    # describe it
    insideInnerCyl.description = 'inside surface of inner quartz cylinder'
    # name its geometry type (options are cylinder, sphere, plane, torus, quadsurface)
    insideInnerCyl.shape = "cylinder"
    # define the shape's parameters
    insideInnerCyl.param_list = [[0, 0, 0], [0, 0, 1], r1[3]]
    # define the inbounds function of the geometry
    insideInnerCyl.inbounds_function = lambda p: np.reshape(
        (p[:, 2, :] < gs.ijar_elevation) and (p[:, 2, :] >= (gs.ijar_elevation - gs.ijar_cyllength)),
        (np.size(p, 0), -1))
    # provide indexes of refraction
    insideInnerCyl.n_outside = gs.n_jar
    insideInnerCyl.n_inside = gs.n_hydraulic
    # provide surface_type and absorption parameters
    insideInnerCyl.surface_type = 'normal'
    insideInnerCyl.absorption = 0
    # add the geometry to the list
    surface_list.append(insideInnerCyl)

    # Define surface #2
    outsideInnerCyl = surface.surface()
    outsideInnerCyl.description = 'outside surface of inner quartz cylinder'
    outsideInnerCyl.shape = "cylinder"
    outsideInnerCyl.param_list = [[0, 0, 0], [0, 0, 1], r1[2]]
    outsideInnerCyl.inbounds_function = lambda p: np.reshape(
        (p[:, 2, :] < gs.ijar_elevation) and (p[:, 2, :] >= (gs.ijar_elevation - gs.ijar_cyllength)),
        (np.size(p, 0), -1))
    outsideInnerCyl.n_outside = gs.n_target
    outsideInnerCyl.n_inside = gs.n_jar
    outsideInnerCyl.surface_type = 'normal'
    outsideInnerCyl.absorption = 0
    surface_list.append(outsideInnerCyl)

    # handle the bubble
    if gs.bubble_present:
        bubble = surface.surface()
        bubble.description = 'bubble'
        bubble.shape = 'sphere'
        bubble.param_list = [gs.bubble_position, gs.bubble_radius]
        bubble.inbounds_function = lambda p: np.reshape(p[:, 2, :] > -500, (np.size(p, 0), -1))
        bubble.n_outside = gs.n_target
        bubble.n_inside = gs.n_air
        bubble.surface_type = 'normal'
        bubble.absorption = 0
        surface_list.append(bubble)

    # generate some light rays 0 this function won't work, it's just a placeholer
    # CreateSomeLightRays DNE
    [raydirections, pixelmap] = CreateSomeLightRays(gs.cam_resolution, cam_pixelpitch, .5 * (1 + gs.cam_resolution),    #What is 'cam_pixelpitch' for? All of these inputs?
                                                    gs.cam_f, gs.cam_pitch + gs.vp_theta - (np.pi / 2), gs.cam_yaw,
                                                    gs.cam_roll, gs.cam_barreld, gs.cam_lenstype)
    # add them to rays container in proper format that will help with calculation
    rays.append(raydirections)
    # do some math and get the ray starting points, and store them in the variable
    ray_startingpoints.append(
        np.matlib.repmat(vp_center + [gs.cam_x, 0, 0] + gs.cam_z * vp_axis + gs.cam_y * np.cross(vp_axis, [1, 0, 0]),
                         (np.size(raydirections, 0), 1)))

    # Return the values
    return [surface_list, rays, ray_startingpoints, pixels]
