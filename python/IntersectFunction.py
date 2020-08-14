import RayToShape

def IntersectFunction(surface, sp, indir):
	#Get the shape and param_list of the surface
	shape = surface.shape
	param_list = surface.param_list

	#Call RayToShape with the correct parameters, and return output
	return RayToShape.RayToShape(shape, sp, indir, param_list)


