# RayTracer2 Python port testing suite
import numpy as np
import math
import random
import RayTracer2
import surface
import matplotlib.pyplot as plt
import sys


#def TestGeometry():
#
#    
#    return ray_startpoints, test_rays, surface_list
    
def main():
        
        
    #test geometry only consists of 1 surface
    surface_list = []
    #sphere
    inner_sphere = surface.surface
    inner_sphere.description = 'sphere centered at origin, radius 10 cm'
    inner_sphere.shape = 'sphere'
    inner_sphere.param_list = [np.array([0,0,0]),10]
    inner_sphere.inbounds_function = lambda p: np.reshape((p[:,2,:] >=-10) * (p[:,2,:]<=10), (p.shape[0],-1))
    inner_sphere.n_outside = np.inf
    inner_sphere.n_inside = 2
    inner_sphere.surface_type = 'normal'
    inner_sphere.absorption = 1
    surface_list.append(inner_sphere)
    
    x=0
    y=0
    z=0
    
    n=10 #starting number of rays
    
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
        
    epsilon = sys.float_info.epsilon
    print("epsilon: " + str(epsilon))

#    ri_data = [] # collect data from each trial
    absorption_data = []
    absorbed_bot = []
    absorbed_top = []



    [ray_interfaces, absorption_table, raytable] = RayTracer2.RayTracer2(ray_startpoints, test_rays, surface_list)
#
#    ri_data.append(ray_interfaces)
    absorption_data.append(absorption_table)

    absorbed_bot.append(sum((np.count_nonzero(x.surface_index == 4) + np.count_nonzero(x.surface_index == (-4)) for x in ray_interfaces))) # num of rays absorbed on bottom
    absorbed_top.append(sum((np.count_nonzero(x.surface_index == 2) + np.count_nonzero(x.surface_index == (-2)) for x in ray_interfaces)))
    report(ray_interfaces, absorption_table, raytable, epsilon)


#    ri_data = np.array(ri_data)
    
    
def report(ray_interfaces, absorption_table, raytable, epsilon): # print relevant info
    for i in range(len(ray_interfaces)):
        print("Scatter # " + str(i+1) + ", # of rays " + str(len(ray_interfaces[i].intersection_point)))
        print("# of times each surface is hit:")
        print("Inner Sphere" + str(np.count_nonzero(ray_interfaces[i].surface_index == 0)))
        #print(ray_interfaces[i].surface_index)

        print("Points of intersection:")
        print(ray_interfaces[i].intersection_point)

        print("Total intensity absorbed by each surface:")
        print("Inner Sphere: " + str(absorption_table[i, 0, 0, :]))
        
        print("Rays escaping geometry: " + str(absorption_table[i,2]))

        print("\n")

if __name__ == "__main__":
    main()

    