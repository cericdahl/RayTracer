#Called by IntersectFunction (defined in IntersectFunction.py) with the proper inputs
def RayToShape(shape, sp, indir, param_list):
    #maybe turn each of these into try-except statements, in case the specified
    # param_list isn't the right size of elements
    if shape == "cylinder":
        output = RayToCylinder(sp, indir, param_list[0], param_list[1], param_list[2])
    elif shape == "sphere":
        output = RayToSphere(sp, indir, param_list[0], param_list[1])
    elif shape == "torus":
        output = RayToTorus(sp, indir, param_list[0], param_list[1], param_list[2], param_list[3])
    elif shape == "plane":
        output = RayToPlane(sp, indir, param_list[0], param_list[1])
    elif shape == "quadsurface":
        output = RayToQuadSurface(sp, indir, param_list[0], param_list[1], param_list[2])
    else:
        raise Exception('Geometry has unrecognized shape')
        
    return output
