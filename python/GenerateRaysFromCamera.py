# % function [rays pixels] = GenerateRaysFromCamera(pinhole_position, resolution, ...
# %     pixel_pitch, pixel_center, focal_length, pitch, yaw, roll, radial_distortion)
# %
# % 12/23/09, CED


# function [ray_direction pixels] = GenerateRaysFromCamera(resolution, pixel_pitch, ...
#     pixel_center, focal_length, pitch, yaw, roll, radial_distortion, lens_type)

def GenerateRaysFromCamera(resolution, pixel_pitch,pixel_center, focal_length, pitch=0, yaw=0, roll=0, radial_distortion=[], lens_type='tan'):

    ##### MAYBE HAVE TO TURN ALL INPUT ARRAYS into numpy arrays to be able to do calculations on it

    #get number of locals
    nargin = len(locals())

    if (nargin<4 or len(focal_length)!=1 or len(yaw)!=1 or len(pitch)!=1 or len(roll)!=1 or len(pixel_center)!=2 or (not bool(pixel_pitch)) or len(pixel_pitch)>2 or len(resolution)!=2):
        raise ValueError('impropper input to GenerateRaysFromCamera')

    if (len(pixel_pitch)==1):
        pixel_pitch = [0, 0] + pixel_pitch

#    %% generate rays in camera-frame
#    % in camera frame, forward is +y, +x is +i, and +z is -j
#    i_pix = repmat((1:resolution(1))',1,resolution(2));
#    j_pix = repmat((1:resolution(2)),resolution(1),1);
    i_pix = np.matlib.repmat(np.transpose(range(resolution[0])),1,resolution[1])
    j_pix = np.matlib.repmat((range(resolution[1])),resolution[0],1)
    pixels = [i_pix[:], j_pix[:]]

    pixel_location = [(pixel_center[0]-i_pix[:])*pixel_pitch[0], np.zeros(len(i_pix),1)-focal_length, -(pixel_center[1]-j_pix[:])*pixel_pitch[1]]

    pixel_d2 = np.sum(np.array(pixel_location[:,[1, 3]])**2,1)

    effective_f = focal_length * (1 + np.sum(np.matlib.repmat(np.transpose(radial_distortion[:]),len(pixel_d2),1)*(np.matlib.repmat(focal_length**(-2) * pixel_d2,1,len(radial_distortion))**np.matlib.repmat(range(len(radial_distortion)),len(pixel_d2),1)), 1))

#    switch lens_type
#        case 'theta'
#            theta = sqrt(pixel_d2)./effective_f;
#        case 'sin'
#            theta = asin(sqrt(pixel_d2)./effective_f);
#        case 'tan'
#            theta = atan(sqrt(pixel_d2)./effective_f);
#        otherwise
#            theta = atan(sqrt(pixel_d2)./effective_f);
#    end

    if (lens_type == 'theta'):
        theta = np.sqrt(pixel_d2)/effective_f
    elif (lens_type ==  'sin'):
        theta = np.arcsin(np.sqrt(pixel_d2)/effective_f)
    elif (lens_type == 'tan'):
        theta = np.arctan(np.sqrt(pixel_d2)/effective_f)
    else:
        theta = np.arctan(np.sqrt(pixel_d2)/effective_f)

    phi = np.arctan2(-pixel_location[:,2],-pixel_location[:,0])
    ray_direction = [np.sin(theta)*np.cos(phi), np.cos(theta), np.sin(theta)*np.sin(phi)]

#    % ray_direction = ray_direction ./ repmat(abs(sqrt(sum(ray_direction.^2,2))),1,3);

#    %% rotate camera
    M1 = [[np.cos(yaw), -np.sin(yaw), 0],[sin(yaw), cos(yaw), 0] [0, 0, 1]]
    M2 = [[1, 0, 0][0, np.cos(pitch), -np.sin(pitch) ][0, np.sin(pitch), np.cos(pitch)]]
    M3 = [[np.cos(roll), 0, np.sin(roll)][ 0, 1, 0 ][-np.sin(roll), 0, np.cos(roll)]]

    M = M1*M2*M3

    ray_direction = np.transpose((M * (np.transpose(ray_direction))))

    return [ray_direction, pixels]
