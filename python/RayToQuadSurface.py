# % RayToQuadsurface specifics:
# % M = 2, rays that miss the sphere give imaginary or NaN outputs for
# % distance_traveled
# %
# % other inputs:
# %       Q                       -   3x3 tensor
# %       P                       -   3 element vector
# %       R                       -   scalar
# %
# % The surface is that which satisfies the equation
# %      x' * Q * x   +   P' * x   +   R   =   0
# %
# % (where here x is 3x1, Q is 3x3, and P is 3x1)
# %
# % the outward-pointing normal at x is taken in the direction of
# %
# %      2 * Q * x   +   P
# %
# %
# % 12/15/09, CED
#
# function [intersection_points surface_normals distance_traveled crossing_into] = ...
#     RayToQuadsurface(starting_points, incoming_directions, Q, P, R)

import sys
import numpy as np
import numpy.matlib

def RayToQuadSurface(starting_points, indir, q, p, r):
    params = len(locals())

    intersection_points = []
    surface_normals = []
    distance_traveled = []
    crossing_into = []

    # check inputs
    if params != 5 or q.size != 9 or len(p) != 3 or r.size != 1 or starting_points.shape[1] != 3 or \
            indir.shape[1] != 3 or len(starting_points) != len(indir) or len(starting_points.shape) != 2 or \
            len(indir.shape) != 2:
        raise Exception('Improper input to RayToQuadsurface')
    q = np.reshape(q, (3, 3))
    p = p[:]
    numrays = starting_points.shape[0]

    """
    # normalize directions
    # incoming_directions(goodray_cut,:) = incoming_directions(goodray_cut,:) ./ ...
    #         repmat(abs(sqrt(sum(incoming_directions.^2,2))),1,3);
    goodray_cut = np.sum(indir ** 2, 1) > 0
    if np.any(goodray_cut):
        indir[goodray_cut, :] = indir[goodray_cut, :] / np.matlib.repmat(np.abs(np.sqrt(np.sum(indir ** 2, 1)))[:, np.newaxis], 1, 3)
    """

    # solve quadratic for distance traveled
    # a*l^2 + b*l + c = 0
    # a = sum((incoming_directions * Q) .* incoming_directions, 2);
    # b = incoming_directions * P + ...
    #     sum((incoming_directions * Q) .* starting_points, 2) + ...
    #     sum((starting_points * Q) .* incoming_directions, 2);
    # c = R + (starting_points * P) + ...
    #     sum((starting_points * Q) .* starting_points, 2);
    a = np.sum((indir * q) * indir, 1)
    b = (indir * p) + np.sum((indir * q) * starting_points, 1) + np.sum((starting_points * q) * indir, 1)
    c = r + (starting_points * p) + np.sum((starting_points * q) * starting_points, 1) # will likely have to turn r into an np array

    # linear cut = a==0 & b!=0
    # the previous line is strictly true, but we also want to avoid rounding error here...
    # linear_cut = b~=0;
    # linear_cut(linear_cut) = abs(4.*a(linear_cut).*c(linear_cut)./(b(linear_cut).*b(linear_cut)))<(100*eps);
    # quad_cut = a~=0 & ~linear_cut;
    # distance_traveled = zeros(numrays,2);
    # distance_traveled(~(linear_cut | quad_cut),:) = NaN;
    linear_cut = np.not_equal(b, 0)
    linear_cut[linear_cut] = np.less(np.abs(4 * a[linear_cut] * c[linear_cut] / (b[linear_cut] ** 2)), (100 * np.spacing(np.float64(1)))) # machine epsilon
    quad_cut = np.logical_and(np.not_equal(a, 0), ~linear_cut)
    distance_traveled = np.zeros(numrays, 2)
    distance_traveled[~(np.logical_or(linear_cut, quad_cut)), :] = np.nan

    # if any(linear_cut)
    #     distance_traveled(linear_cut,1) = -c(linear_cut) ./ b(linear_cut);
    #     distance_traveled(linear_cut,2) = -b(linear_cut) ./ a(linear_cut);
    # if any(quad_cut)
    #     distance_traveled(quad_cut,:) = ...
    #         repmat(-.5 .* b(quad_cut) ./ a(quad_cut),1,2) + (.5 .* sqrt(b(quad_cut).^2 - 4 .* a(quad_cut) .* c(quad_cut)) ./ a(quad_cut)) * [1 -1];
    if np.any(linear_cut):
        distance_traveled[linear_cut, 0] = -c[linear_cut] / b[linear_cut]
        distance_traveled[linear_cut, 1] = -b[linear_cut] / a[linear_cut]
    if np.any(quad_cut):
        distance_traveled[quad_cut, :] = (np.matlib.repmat(-0.5 * b[quad_cut] / a[quad_cut], 1, 2) +
            (0.5 * np.sqrt(b[quad_cut]**2 - 4 * a[quad_cut] * c[quad_cut]) / a[quad_cut])) * [1, -1] # check repmat shape

    # find intersection points
    intersection_points = starting_points[:,:,np.newaxis] + distance_traveled[:, np.newaxis, :] * indir[:, :, np.newaxis] # No need to reshape or tile starting_points or indir, will be handled automatically by numpy broadcasting

    # find surface normals
    # surface_normals = zeros(size(intersection_points));
    # for n=1:size(surface_normals,3)
    #     surface_normals(:,:,n) = 2 * intersection_points(:,:,n) * Q + repmat(P', numrays, 1);
    #     goodnormal_cut = sum(surface_normals(:,:,n).^2, 2)>0;
    #     surface_normals(goodnormal_cut,:,n) = surface_normals(goodnormal_cut,:,n) ./ ...
    #         repmat(abs(sqrt(sum(surface_normals(goodnormal_cut,:,n).^2, 2))),1,3);
    surface_normals = np.zeros(intersection_points.shape)
    for n in range(0, surface_normals.shape[2] + 1):
        surface_normals[:, :, n] = 2 * intersection_points[:, :, n] * q + np.matlib.repmat(np.transpose(p), numrays, 1) # check repmat shape
        goodnormal_cut = np.sum(surface_normals[:, :, n] ** 2, 1) > 0
        surface_normals[goodnormal_cut, :, n] = surface_normals[goodnormal_cut, :, n] / \
            np.matlib.repmat(np.abs(np.sqrt(np.sum(surface_normals[goodnormal_cut, :, n] ** 2, 1))), 1, 3) # check repmat shape

    # crossing_into = round(-sign(sum(repmat(incoming_directions,[1,1,2]) .* surface_normals,2)));
    # surface_normals = surface_normals .* repmat(crossing_into,[1 3 1]);
    # crossing_into = reshape(crossing_into,[],2);
    crossing_into = np.round_(np.sign(np.sum(indir[:, :, np.newaxis] * surface_normals, axis=1)))
    surface_normals = surface_normals * crossing_into[:, np.newaxis, :]
    crossing_into = np.reshape(crossing_into, (-1, 2))

    return [intersection_points, surface_normals, distance_traveled, crossing_into]