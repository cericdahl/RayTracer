# % function [intersection_points surface_normals distance_traveled crossing_into] = ...
# %     RayToCylinder(starting_points, incoming_directions, cylinder_center, cylinder_axis, cylinder_radius)
# % 
# % One of several RayToXXXXX functions used for RayTracing
# %
# % All RayToXXXXX functions have the following inputs and outputs in common:
# %
# %   inputs:
# %       starting_points         -   N-by-3 array, where N is the number of
# %                                     rays being followed.  Gives the
# %                                     starting position for each ray
# %       incoming_directions     -   N-by-3 array, gives the
# %                                     forward-direction of each ray
# %   
# %   outputs:
# %       intersection_points     -   N-by-3-by-M array of intersection
# %                                     points, M varies between functions
# %       surface_normals         -   N-by-3-by-M array of surface normals.
# %                                     Normals point opposite to
# %                                     incoming_directions (dot product < 0)
# %       distance_traveled       -   N-by-M array of distances (negative for
# %                                     hits behind the starting point, NaN,
# %                                     +inf, -inf, or imaginary values are
# %                                     non-intersections)
# %       crossing_into           -   N-by-M array of +1,0,-1, where +1
# %                                     indicates the ray crossing into the
# %                                     surface, (surface_normal is aligned
# %                                     with outward-normal), -1 indicates
# %                                     the ray is leaving the surface
# %                                     (surface_normal is anti-aligned with
# %                                     outward-normal), and 0 indicates a
# %                                     glancing blow
# %
# % RayToCylinder specifics:
# % M = 2, rays that miss the cylinder give imaginary outputs for
# % distance_traveled
# %
# % other inputs:
# %       cylinder_center         -   3-element vector, giving a point along
# %                                     the cylinder axis
# %       cylinder_axis           -   3-element vector pointing along the
# %                                     cylinder axis
# %       cylinder_radius         -   scalar giving cylinder radius
# %                                    
# %
# % 12/15/09, CED

# function [intersection_points surface_normals distance_traveled crossing_into] = ...
    # RayToCylinder(starting_points, incoming_directions, cylinder_center, cylinder_axis, cylinder_radius)

import numpy as np

def RayToCylinder(starting_points, indir, cyl_center, cyl_axis, cylinder_radius):

    intersection_points = []
    surface_normals = []
    distance_traveled = []
    crossing_into = []

    # %% check inputs
    if len(cylinder_center)!=3 or len(cylinder_axis)!=3 or len(cylinder_radius)!=1 or np.size(starting_points,1)!=3 or np.size(incoming_directions,1)!=3 or np.size(starting_points,0)!=np.size(incoming_directions,0):
        raise Exception('Impropper input to RayToCylinder')
    cylinder_center = np.transpose(cyl_center)
    cylinder_axis = np.transpose(cyl_axis)
    numrays = np.size(starting_points,0)

    # %% normalize directions
    goodray_cut = np.sum(incoming_directions**2,1)>0
    if any(goodray_cut):
        incoming_directions[goodray_cut,:] = indir[goodray_cut,:] / np.matlib.repmat(abs(np.sqrt(np.sum(incoming_directions**2,1))),1,3)
    if np.sum(cylinder_axis**2) > 0:
        cylinder_axis = cylinder_axis / abs(np.sqrt(np.sum(cylinder_axis**2)))
    else:
        raise Exception('Invalid cylinder for RayToCylinder')

    # %% solve quadratic for distance_traveled
    # % || u + l*v || = r
    # % a*l^2 + b*l + c = 0, a = 1

    x = np.array(starting_points) - np.matlib.repmat(cylinder_center,numrays,1)
    u = (x * np.transpose(cylinder_axis) * cylinder_axis) - x


    v = (np.array(incoming_directions) * np.transpose(cylinder_axis) * cylinder_axis) - np.array(incoming_directions)

    a = np.sum(np.multiply(v,v), 1)
    b = 2 * np.sum(np.multiply(u,v), 1)
    c = np.sum(np.multiply(u,u), 1) - cylinder_radius**2

    linear_cut = (a==0 and b!=0)
    quad_cut = (a!=0)
    distance_traveled = np.zeros((numrays,2))


    ############# BOOKMARK ######## 
    # distance_traveled(~(linear_cut | quad_cut),:) = NaN;
    cutIndex = np.logical_or(linear_cut, quad_cut)
    distance_traveled[np.logical_not(cutIndex)] = np.nan

    if any(linear_cut): #Might have to switch dimensions linear_cut and :
        distance_traveled[linear_cut,:] = np.matlib.repmat(-c[linear_cut] / b[linear_cut], 1, 2)
    if any(quad_cut):
        distance_traveled[quad_cut,:] = np.matlib.repmat(-.5 * b[quad_cut] / a[quad_cut],1,2) + (.5 * np.sqrt(b[quad_cut]**2 - 4 * a[quad_cut] * c[quad_cut]) / a[quad_cut]) * np.array([1, -1])

    # %% find intersection_points
    # intersection_points = repmat(starting_points,[1,1,2]) + repmat(reshape(distance_traveled,[],1,2),[1,3,1]) .* repmat(incoming_directions,[1,1,2]);
    #The first reshape might be wrong, not sure how to reshape into 3d array. Should be right according to this link: https://bic-berkeley.github.io/psych-214-fall-2016/reshape_and_3d.html
    intersection_points = np.tile(starting_points,[1,1,2]) + np.tile(np.reshape(distance_traveled,(-1,1,2)),[1,3,1]) * np.matlib.repmat(incoming_directions,[1,1,2])


    # %% find surface_normals
    x = intersection_points[:,:,0] - np.matlib.repmat(cylinder_center,(numrays,1))
    u = (x * np.transpose(cylinder_axis) * cylinder_axis) - x
    y = intersection_points[:,:,1] - np.matlib.repmat(cylinder_center,(numrays,1))
    v = (y * np.transpose(cylinder_axis) * cylinder_axis) - y
    surface_normals = np.zeros(np.size(intersection_points))
    surface_normals[:,:,1] = u / cylinder_radius
    surface_normals[:,:,2] = v / cylinder_radius
    #Not sure what the [] is at the end of repmat
    # crossing_into = round(sign(sum(repmat(incoming_directions,[1,1,2]) .* surface_normals,2)));
    # surface_normals = -surface_normals .* repmat(crossing_into,[1 3 1]);
    crossing_into = np.round(np.sign(np.sum(np.matlib.repmat(incoming_directions,1,1,2) * surface_normals,2)))
    surface_normals = -surface_normals * np.matlib.repmat(crossing_into,1, 3, 1)

    crossing_into = np.reshape(crossing_into,(-1,2))

    #Return a list of the return values (python)
    return [intersection_points, surface_normals, distance_traveled, crossing_into]



