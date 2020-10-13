# RayTracer2 Python port testing suite
import numpy as np
import math
import random
import RayTracer2
import surface
import matplotlib.pyplot as plt
import sys


def TestGeometry(z_in):
    # First, create ray starting points and isotropic rays
    # Coordinates of startingpoint -- same for all rays
    x = 0
    y = 0
    z = z_in

    n = 10 # number of rays

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
        polar = np.arccos(a) # -1 < theta < 1 so random numbers are not biased towards the pole

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
    bot_cyl.inbounds_function = lambda p: np.reshape((p[:, 2, :] > 0) * (p[:, 2, :] < 5), (p.shape[0], -1))
    bot_cyl.n_outside = 1.5
    bot_cyl.n_inside = 1.5
    bot_cyl.surface_type = 'unified'
    bot_cyl.unifiedparams = [0, 0, 0, 0, 0]
    bot_cyl.absorption = 0
    surface_list.append(bot_cyl)

    # Top cylinder
    top_cyl = surface.surface()
    top_cyl.description = 'top cylinder along z-axis, 10cm radius from z=5 to z=10'
    top_cyl.shape = 'cylinder'
    top_cyl.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), 10]
    top_cyl.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= 5) * (p[:, 2, :] < 10), (p.shape[0], -1))
    top_cyl.n_outside = 1.5
    top_cyl.n_inside = 2
    top_cyl.surface_type = 'unified'
    top_cyl.unifiedparams = [0, 0, 0, 0, 0]
    top_cyl.absorption = 0
    surface_list.append(top_cyl)

    # Top cap
    top = surface.surface()
    top.description = 'top cap, disk centered on z-axis with radius 10 and z=10'
    top.shape = 'plane'
    top.param_list = [np.array([0, 0, 10]), np.array([0, 0, 1])]
    top.inbounds_function = lambda p: np.reshape((p[:, 0] ** 2 + p[:, 1] ** 2) < 100, (p.shape[0], -1))
    # Direction of normal vector considered 'outside'
    top.n_outside = 1.5
    top.n_inside = 2
    top.surface_type = 'normal'
    top.absorption = 1
    surface_list.append(top)

    # Middle
    mid = surface.surface()
    mid.description = 'middle disk centered on z-axis with radius 10 and z=5'
    mid.shape = 'plane'
    mid.param_list = [np.array([0, 0, 5]), np.array([0, 0, 1])]
    mid.inbounds_function = lambda p: np.reshape((p[:, 0] ** 2 + p[:, 1] ** 2 < 100), (p.shape[0], -1))
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
    bottom.inbounds_function = lambda p: np.reshape((p[:, 0] ** 2 + p[:, 1] ** 2 < 100), (p.shape[0], -1))
    bottom.n_outside = 1.5
    bottom.n_inside = 1.5
    bottom.surface_type = 'normal'
    bottom.absorption = 1
    surface_list.append(bottom)

    return ray_startpoints, test_rays, surface_list


def main():
    epsilon = sys.float_info.epsilon
    print("epsilon: " + str(epsilon))

    ri_data = [] # collect data from each trial
    absorption_data = []
    absorbed_bot = []
    absorbed_top = []

    """
    for i in range(1): # test a bunch of times
        for z in np.arange(.2, .4, .2): # move z up the center of cylinder in steps of 0.2
            [starts, rays, surfaces] = TestGeometry(z)

            [ray_interfaces, absorption_table, raytable] = RayTracer2.RayTracer2(starts, rays, surfaces)

            ri_data.append(ray_interfaces)
            absorption_data.append(absorption_table)

            absorbed_bot.append(sum((np.count_nonzero(x.surface_index == 4) + np.count_nonzero(x.surface_index == (-4)) for x in ray_interfaces))) # num of rays absorbed on bottom
            absorbed_top.append(sum((np.count_nonzero(x.surface_index == 2) + np.count_nonzero(x.surface_index == (-2)) for x in ray_interfaces)))

            print("z = " + str(z))
            report(ray_interfaces, absorption_table, raytable, epsilon)


    ri_data = np.array(ri_data)
    """

    [starts, rays, surfaces] = TestGeometry(7)
    [ray_interfaces, absorption_table, raytable] = RayTracer2.RayTracer2(starts, rays, surfaces)
    report(ray_interfaces, absorption_table, raytable, epsilon)


    # # plot
    # x_data = np.arange(.2, 10, .2)
    # width = 0.05
    # labels = np.arange(.2, 10, .2).round(decimals=1) # fix x-labels
    #
    # fig, ax = plt.subplots()
    #
    # ax.bar(x_data, absorbed_bot, width, label='absorbed bottom')
    # ax.bar(x_data, absorbed_top, width, bottom=absorbed_bot, label='absorbed top')
    # plt.xticks(np.arange(.2, 10, .2), labels, size='small', rotation='vertical')
    #
    # ax.set_ylabel('# of rays absorbed')
    # ax.set_xlabel('z-start')
    # ax.set_title('# of rays absorbed at caps by starting height')
    # ax.legend()
    #
    # plt.show()




def report(ray_interfaces, absorption_table, raytable, epsilon): # print relevant info
    for i in range(len(ray_interfaces)):
        print("Scatter # " + str(i+1) + ", # of rays " + str(len(ray_interfaces[i].intersection_point)))
        print("# of times each surface is hit:")
        print("Bot Cyl: " + str(np.count_nonzero(ray_interfaces[i].surface_index == 0)))
        print("Top Cyl: " + str(np.count_nonzero(ray_interfaces[i].surface_index == 1) + np.count_nonzero(ray_interfaces[i].surface_index == (-1))))
        print("Top Cap: " + str(np.count_nonzero(ray_interfaces[i].surface_index == 2) + np.count_nonzero(ray_interfaces[i].surface_index == (-2))))
        print("Mid Interface: " + str(np.count_nonzero(ray_interfaces[i].surface_index == 3) + np.count_nonzero(ray_interfaces[i].surface_index == (-3))))
        print("Bot Cap: " + str(np.count_nonzero(ray_interfaces[i].surface_index == 4) + np.count_nonzero(ray_interfaces[i].surface_index == (-4))))
        #print(ray_interfaces[i].surface_index)

        print("Points of intersection:")
        print(ray_interfaces[i].intersection_point)

        print("Total intensity absorbed by each surface:")
        print("Bot Cyl: " + str(absorption_table[i, 0, 0, :]))
        #if absorption_table[i, 0, 0, 0] > (5 * epsilon) or absorption_table[i, 0, 0, 0] < -(5 * epsilon) or absorption_table[i, 0, 0, 1] > (5 * epsilon) or absorption_table[i, 0, 0, 1] < -(5 * epsilon):
            #raise Exception("Bottom Cylinder absorbed!!!")
        print("Top Cyl: " + str(absorption_table[i, 0, 1, :]))
        print("Top Cap: " + str(absorption_table[i, 0, 2, :]))
        print("Mid Interface: " + str(absorption_table[i, 0, 3, :]))
        print("Bot Cap: " + str(absorption_table[i, 0, 4, :]))

        print("Rays escaping geometry: " + str(absorption_table[i,2]))

        print("\n")


if __name__ == "__main__":
    main()

