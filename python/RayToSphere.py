# % RayToSphere specifics:
# % M = 2, rays that miss the sphere give imaginary outputs for
# % distance_traveled
# %
# % other inputs:
# %       sphere_center           -   3-element vector, giving the center of
# %                                     the sphere
# %       sphere_radius           -   scalar giving sphere radius
# %
# %
# % 12/15/09, CED
# function [intersection_points surface_normals distance_traveled crossing_into] = ...
#     RayToSphere(starting_points, incoming_directions, sphere_center, sphere_radius)

import sys
import numpy as np

def RayToSphere(starting_points, indir, sphere_center, sphere_radius):
    params = len(locals())

    #print("sphere!!")

    intersection_points = []
    surface_normals = []
    distance_traveled = []
    crossing_into = []

    # check inputs
    if params != 4 or len(sphere_center) != 3 or np.size(sphere_radius) != 1 or starting_points.shape[1] != 3 or \
            indir.shape[1] != 3 or starting_points.shape[0] != indir.shape[0]:
        raise Exception('Improper input to RayToSphere')
    sphere_center = np.transpose(sphere_center[:,np.newaxis])
    numrays = starting_points.shape[0]

    """
    # normalize directions
    goodray_cut = np.sum(indir ** 2, 1) > 0
    if np.any(goodray_cut):
        indir[goodray_cut, :] = indir[goodray_cut, :] / np.matlib.repmat(np.abs(np.sqrt(np.sum(indir ** 2, 1)))[:, np.newaxis], 1, 3)
    """

    # solve quadratic for distance traveled
    a = np.sum(indir ** 2, axis=1)
    b = 2 * np.sum(indir * (starting_points - np.matlib.repmat(sphere_center, numrays, 1)), axis=1)
    c = np.sum((starting_points - np.matlib.repmat(sphere_center, numrays, 1)) ** 2, axis=1) - (sphere_radius ** 2)

    distance_traveled = (np.matlib.repmat((-0.5 * b / a)[:,np.newaxis], 1, 2) + (0.5 * np.sqrt(b ** 2 - 4 * a * c) / a)[:,np.newaxis]) * np.array([1, -1]) # correct
    #print("distance: " + str(distance_traveled))

    # find intersection points
    intersection_points = starting_points[:, :, np.newaxis] + distance_traveled[:, np.newaxis, :] * indir[:, :, np.newaxis]
    #print("intersection: " + str(intersection_points))
    #print(intersection_points.shape)

    # find surface normals
    # surface_normals = (intersection_points - repmat(sphere_center,[numrays,1,2])) ./ sphere_radius;
    # crossing_into = round(-sign(sum(repmat(incoming_directions,[1,1,2]) .* surface_normals,2)));
    # surface_normals = surface_normals .* repmat(crossing_into,[1 3 1]);
    # crossing_into = reshape(crossing_into,[],2);
    #surface_normals = (intersection_points - np.tile(sphere_center[:, np.newaxis], (numrays, 1, 2))) / sphere_radius # np.tile or [:, :, np.newaxis]?
    surface_normals = (intersection_points - sphere_center[:, :, np.newaxis]) / sphere_radius
    crossing_into = np.round_(-np.sign(np.sum(indir[:, :, np.newaxis] * surface_normals, axis=1)))
    surface_normals = surface_normals * crossing_into[:, np.newaxis, :]
    crossing_into = np.reshape(crossing_into, (-1, 2))

    return [intersection_points, surface_normals, distance_traveled, crossing_into]

