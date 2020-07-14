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

class surface(object):
    def __init__(self):
        # attributes
        self.description = None
        self.n_outside = None
        self.n_inside = None
        self.surface_type = None
        self.absorption = None
        # Add these to replace assigning RayTo____() function at runtime
        self.shape = None  # For now they will be strings, later maybe change to integers to represent each type
        self.param_list = None
        self.inbounds_function = None

        # No longer have these
        # self.intersect_function= None
