import numpy as np
import math
import RayTracer2
import surface
import matplotlib.pyplot as plt

# Check retroreflection in UnifiedReflectorModel with a retroreflecting panel covered in a hemisphere (for collecting
# distribution). Fire light straight down onto the center of the panel and observe distribution on hemisphere.


def main():
    surface_list = []

    # hemisphere
    hemi = surface.surface()
    hemi.description = '10cm radius hemisphere, in positive z'
    hemi.shape = 'sphere'
    hemi.param_list = [np.array([0, 0, 0]), 10]
    hemi.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= 0) * (p[:, 2, :] <= 10), (p.shape[0], -1))
    hemi.n_outside = np.inf
    hemi.n_inside = 1.5
    hemi.surface_type = 'normal'
    hemi.absorption = 1
    surface_list.append(hemi)

    # panel
    panel = surface.surface()
    panel.description = 'retroreflecting panel cut to 10cm radius disk on xy-plane'
    panel.shape = 'plane'
    panel.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1])]
    panel.inbounds_function = lambda p: np.reshape((p[:, 0] ** 2 + p[:, 1] ** 2) < 100, (p.shape[0], -1))
    panel.n_outside = 1.5
    panel.n_inside = np.inf
    panel.surface_type = 'unified'
    panel.unifiedparams = [0, 0, .25, .5, .1]
    panel.absorption = 0
    surface_list.append(panel)
    print("size test: " + str(np.size(panel.unifiedparams)))

    # construct light rays
    # start them along z-axis, pointing down for normal incidence
    # start at (-5, 0, 5), pointing towards origin for 45 degree incidence
    x = -5
    y = 0
    z = 5

    n = 1000000  # number of rays

    ray_startpoints = np.empty((n, 3))
    ray_startpoints[..., 0] = x
    ray_startpoints[..., 1] = y
    ray_startpoints[..., 2] = z

    # rays = initial [forward direction (3), s1 polarization axis (3), s0, s1, s2, s3]
    # Let s0 be 1 and rays be unpolarized (other Stokes parameters = 0 + no polarization axis)
    test_rays = np.zeros((n, 10))
    test_rays[..., 3] = 1
    test_rays[..., 6] = 1

    test_rays[:, 0] = 1
    test_rays[:, 1] = 0
    test_rays[:, 2] = -1

    [ray_interfaces, absorption_table, raytable] = RayTracer2.RayTracer2(ray_startpoints, test_rays, surface_list)


    # PLOT
    # calculate spherical coordinates
    x_int = ray_interfaces[1].intersection_point[:, 0]
    y_int = ray_interfaces[1].intersection_point[:, 1]
    z_int = ray_interfaces[1].intersection_point[:, 2]

    theta = np.arctan2(np.sqrt(x_int**2 + y_int**2), z_int) # 0 < theta < pi/2
    phi = np.arctan2(y_int, x_int) # -pi < phi < pi
    phi[np.isnan(phi)] = 0
    phi = (phi + 2 * math.pi) % (2 * math.pi) # convert to 0 < phi < 2pi

    # spherical coords of initial rays
    theta_0 = np.matlib.repmat(np.cos(np.arctan2(np.sqrt(x ** 2 + y ** 2), z)), len(theta), 1)
    phi_0 = np.matlib.repmat(np.arctan2(y, x), len(theta), 1)
    phi_0[np.isnan(phi_0)] = 0
    #phi_0 = (phi_0 + 2 * math.pi) % (2 * math.pi)

    # Perfect Reflection Check
    points = np.concatenate((phi[:, np.newaxis], np.cos(theta)[:, np.newaxis]), axis=1)
    counts = np.equal(points, np.concatenate((phi_0, np.cos(theta_0)), axis=1)) # cos(theta), phi the same

    fig, ax = plt.subplots()
    ax.scatter(phi, np.cos(theta))

    plt.xlim(0, 2 * math.pi)
    plt.ylim(0, 1)

    ax.set_ylabel("cos\u03B8") # theta
    ax.set_xlabel("\u03C6  (azimuth)") # phi
    ax.grid(True)

    print("intersection: " + str(ray_interfaces[1].intersection_point))
    print("points: " + str(points))
    print("counts: " + str(counts))
    print("check: " + str(np.concatenate((phi_0, np.cos(theta_0)), axis=1)))

    plt.text(.96, .08, "# rays perfectly reflected = {}".format(int(np.sum(counts))), bbox={'facecolor':'w','pad':5}, ha="right", va="top", transform=plt.gca().transAxes) #check the counts/2

    plt.show()



if __name__ == "__main__":
    main()

