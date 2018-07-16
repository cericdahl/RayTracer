# RayTracer
Optical ray tracing code

Propagate rays through a surface based geometry.

Surfaces are defined by their geometry and optical properties.  Supported geometries include sections of linear, quadratic, and toroidal surfaces, with easy implementation shortcuts for planes, cylinders, and spheres.  Supported optical properties include dielectric interfaces, diffuse and retro-reflectors, and a unified absorber/reflector model based on the Geant4 "UNIFIED" optical photon model.

Bulk properties (index of refraction, absorption and Rayleigh scattering) are implemented at the surfaces (i.e. you define the bulk on either side of each surface -- no checks for inconsistent geometries).

Rays are defined by their propagation direction and stokes parameters (intensity and polarization) and are followed until:
 - intensity drops below threshold
 - ray leaves geometry
 - maximum number of scatters is exceeded
 
For geometries without random scattering (i.e. no Rayleigh scattering and no diffuse reflectors) the program can be run with or without Monte Carlo style dice rolling.  E.g., at a dielectric interface, you can choose to follow both the reflected and refracted ray (so the number of tracked rays increases as the simulation proceeds), or to roll the dice according to the relative intensities and follow only one ray.  For diffuse reflection or Rayleigh scattering (or any other process where there is a continuous angular distribution) only Monte Carlo style simulation is possible.

Supported platforms:
 - MATLAB
 - working on a numpy version
