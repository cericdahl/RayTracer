# define ray_interfaces struct
# each rayInterfaces object contains all of the data for one scatter
# each attribute is an nparray
class rayInterfaces(object):
    def __init__(self):
        # attributes
        self.incoming_ray = None
        self.reflected_ray = None
        self.refracted_ray = None
        self.intersection_point = None
        self.surface_normal = None
        self.ray_index = None
        self.surface_index = None
        self.distance_traveled = None
        self.n_incident = None
        self.n_transmitted = None
        self.bulkabs_incident = None
        self.bulkabs_transmitted = None
        self.rayleigh_incident = None
        self.rayleigh_transmitted = None
