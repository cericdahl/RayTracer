# RayTracer2 Python port testing suite
import numpy as np
import math
import random
import RayTracer2
import surface


def TestGeometry():
    # First, create ray starting points and isotropic rays
    # Coordinates of startingpoint -- is the same for all rays
    x = 0
    y = 0
    z = 2

    n = 1000  # number of rays

    ray_startpoints = np.empty((n, 3))
    ray_startpoints[..., 0] = x
    ray_startpoints[..., 1] = y
    ray_startpoints[..., 2] = z

    # rays = initial [forward direction (3), s1 polarization axis (3), s0, s1, s2, s3]
    # Let s0 be 1 and rays be unpolarized (other Stokes parameters = 0 + no polarization axis)
    test_rays = np.zeros((n, 10))
    test_rays[..., 3] = 1
    test_rays[..., 6] = 1

    # Assign initial forward directions in spherical coords for easier isotropism
    for i in range(n):
        azimuth = random.random() * 2 * math.pi
        a = np.random.uniform(-1, 1)
        polar = np.arccos(a) # -1 < theta < 1 so random numbers are not biased towards the poles

        test_rays[i, 0] = np.sin(polar) * np.cos(azimuth) # x
        test_rays[i, 1] = np.sin(polar) * np.sin(azimuth) # y
        test_rays[i, 2] = np.cos(polar)  # z

    # Test geometry consists of 5 surfaces
    surface_list = []
    # Bottom cylinder
    bot_cyl = surface.surface()
    bot_cyl.description = 'bottom cylinder along z-axis, 10cm radius from z=0 to z=5'
    bot_cyl.shape = 'cylinder'
    bot_cyl.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), 10]
    bot_cyl.inbounds_function = lambda p: np.reshape((p[:, 2, :] > 0) * (p[:, 2, :] < 5), (np.size(p, 0), -1)) # Replaced 3 w/ 2
    bot_cyl.n_outside = np.inf
    bot_cyl.n_inside = 1.5
    bot_cyl.surface_type = 'normal'
    bot_cyl.absorption = 0
    surface_list.append(bot_cyl)

    # Top cylinder
    top_cyl = surface.surface()
    top_cyl.description = 'top cylinder along z-axis, 10cm radius from z=5 to z=10'
    top_cyl.shape = 'cylinder'
    top_cyl.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), 10]
    top_cyl.inbounds_function = lambda p: np.reshape((p[:, 2, :] > 5) * (p[:, 2, :] < 10), (np.size(p, 0), -1))
    top_cyl.n_outside = np.inf
    top_cyl.n_inside = 2
    top_cyl.surface_type = 'normal'
    top_cyl.absorption = 0
    surface_list.append(top_cyl)

    # Top cap
    top = surface.surface()
    top.description = 'top cap, disk centered on z-axis with radius 10 and z=10'
    top.shape = 'plane'
    top.param_list = [np.array([0, 0, 10]), np.array([0, 0, 1])]
    top.inbounds_function = lambda p: np.reshape((p[:, 0] ** 2 + p[:, 1] ** 2) < 100, (p.shape[0], -1)) # shifted indices down 1 and eliminated 3rd dimension index
    # Direction of normal vector considered 'outside'
    top.n_outside = np.inf
    top.n_inside = 2
    top.surface_type = 'normal'
    top.absorption = 1
    surface_list.append(top)

    # Middle
    mid = surface.surface()
    mid.description = 'middle, disk centered on z-axis with radius 10 and z=5'
    mid.shape = 'plane'
    mid.param_list = [np.array([0, 0, 5]), np.array([0, 0, 1])]
    mid.inbounds_function = lambda p: np.reshape((p[:, 0] ** 2 + p[:, 1] ** 2 < 100), (np.size(p, 0), -1))
    mid.n_outside = 2
    mid.n_inside = 1.5
    mid.surface_type = 'normal'
    mid.absorption = 0
    surface_list.append(mid)

    # Bottom cap
    bottom = surface.surface()
    bottom.description = 'bottom cap, disk centered on z-axis with radius 10 and z=0'
    bottom.shape = 'plane'
    bottom.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1])]
    bottom.inbounds_function = lambda p: np.reshape((p[:, 0] ** 2 + p[:, 1] ** 2 < 100), (np.size(p, 0), -1))
    bottom.n_outside = 1.5
    bottom.n_inside = np.inf
    bottom.surface_type = 'normal'
    bottom.absorption = 1
    surface_list.append(bottom)

    return ray_startpoints, test_rays, surface_list


def main():
    [starts, rays, surfaces] = TestGeometry()

    """RAYTOCYLINDER BEING CALLED A THIRD TIME AFTER ALL SURFACES -- WHY?"""
    # return of RayTracer2: [ray_interfaces, absorption_table, raytable]

    [ray_interfaces, absorption_table, raytable] = RayTracer2.RayTracer2(starts, rays, surfaces)

    print("ray_interfaces:")
    print(ray_interfaces.shape)
    print("absorption_table:")
    print(absorption_table)
    print("raytable:")
    print(raytable)

if __name__ == "__main__":
    main()

