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
import numpy.matlib

def RayToCylinder(starting_points, indir, cyl_center, cyl_axis, cylinder_radius):

    intersection_points = []
    surface_normals = []
    distance_traveled = []
    crossing_into = []

    # %% check inputs                                        |  Check this       |
    if len(cyl_center) != 3 or len(cyl_axis) != 3 or cylinder_radius == 0 or starting_points.shape[1] != 3 or np.size(indir,1) != 3 or starting_points.shape[0] != indir.shape[0]:
        raise Exception('Improper input to RayToCylinder')
    cylinder_center = np.transpose(cyl_center) # 1x3 --> 3x1 (3,)
    cylinder_axis = np.transpose(cyl_axis)
    numrays = starting_points.shape[0] # numrays = n

    # %% normalize directions
    goodray_cut = np.sum(indir**2, 1) > 0
    if any(goodray_cut):
        indir[goodray_cut,:] = indir[goodray_cut,:] / np.transpose(np.matlib.repmat(abs(np.sqrt(np.sum(indir**2,1))), 3, 1)) # Using goodray_cut for indexing applies normalization to all good rays | was: ValueError: operands could not be broadcast together with shapes (1000,3) (1,3000)
    if np.sum(cylinder_axis**2) > 0:
        cylinder_axis = cylinder_axis / abs(np.sqrt(np.sum(cylinder_axis**2)))
    else:
        raise Exception('Invalid cylinder for RayToCylinder')

    # %% solve quadratic for distance_traveled
    # % || u + l*v || = r
    # % a*l^2 + b*l + c = 0, a = 1

    x = np.array(starting_points) - np.matlib.repmat(cylinder_center,numrays,1)
    u = (x * np.transpose(cylinder_axis) * cylinder_axis) - x


    v = (np.array(indir) * np.transpose(cylinder_axis) * cylinder_axis) - np.array(indir)

    a = np.sum(np.multiply(v,v), 1)
    b = 2 * np.sum(np.multiply(u,v), 1)
    c = np.sum(np.multiply(u,u), 1) - cylinder_radius**2

    linear_cut = (np.logical_and(np.equal(a, 0), np.not_equal(b, 0)))
    quad_cut = (np.not_equal(a, 0))
    distance_traveled = np.zeros((numrays,2))


    ############# BOOKMARK ######## 
    # distance_traveled(~(linear_cut | quad_cut),:) = NaN;
    cutIndex = np.logical_or(linear_cut, quad_cut)
    distance_traveled[np.logical_not(cutIndex)] = np.nan

    if any(linear_cut): #Might have to switch dimensions linear_cut and :
        distance_traveled[linear_cut,:] = np.matlib.repmat(-c[linear_cut] / b[linear_cut], 2, 1) # Swapped 2 and 1
    if any(quad_cut):
        distance_traveled[quad_cut,:] = np.add(np.transpose(np.matlib.repmat(-.5 * b[quad_cut] / a[quad_cut], 2, 1)), (.5 * np.sqrt(b[quad_cut]**2 - 4 * a[quad_cut] * c[quad_cut]) / a[quad_cut])[:, np.newaxis] * np.array([1, -1])) # Added [:, np.newaxis], swapped 2 and 1 in first term, transposed first term, and used np.add instead of +


    # %% find intersection_points
    # intersection_points = repmat(starting_points,[1,1,2]) + repmat(reshape(distance_traveled,[],1,2),[1,3,1]) .* repmat(incoming_directions,[1,1,2]);
    intersection_points = starting_points[:,:,np.newaxis] + distance_traveled[:, np.newaxis, :] * indir[:, :, np.newaxis] # No need to reshape or tile starting_points or indir, will be handled automatically by numpy broadcasting


    # %% find surface_normals
    x = intersection_points[:,:,0] - np.matlib.repmat(cylinder_center, numrays, 1) # intersection_points[:,:,0] takes the first set of quadratic solutions
    u = (x * np.transpose(cylinder_axis) * cylinder_axis) - x
    y = intersection_points[:,:,1] - np.matlib.repmat(cylinder_center, numrays, 1)
    v = (y * np.transpose(cylinder_axis) * cylinder_axis) - y
    surface_normals = np.zeros(np.shape(intersection_points)) # Replaced 'size' w/ 'shape'
    surface_normals[:,:,0] = u / cylinder_radius
    surface_normals[:,:,1] = v / cylinder_radius
    #Not sure what the [] is at the end of repmat
    # crossing_into = round(sign(sum(repmat(incoming_directions,[1,1,2]) .* surface_normals,2)));
    # surface_normals = -surface_normals .* repmat(crossing_into,[1 3 1]);
    crossing_into = np.round_(np.sign(np.sum(indir[:,:,np.newaxis] * surface_normals,axis=1)))
    surface_normals = -surface_normals * crossing_into[:,np.newaxis,:]

    crossing_into = np.reshape(crossing_into,(-1,2))

    #Return a list of the return values (python)
    return [intersection_points, surface_normals, distance_traveled, crossing_into]



