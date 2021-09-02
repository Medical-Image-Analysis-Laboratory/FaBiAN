%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that compiles the name of the output folder where the         %
%  simulated images are saved according to the gestational age of the     %
%  fetus, to the level of fetal movements, to the acquisition plane and   %
%  to the displacement of the fetal brain volume in a given orientation   %
%  from one series to the other.                                          %
%                                                                         %
%                output_folder = output_name(          GA, ...            %
%                                            motion_level, ...            %
%                                             orientation, ...            %
%                                                shift_mm);               %
%                                                                         %
%  inputs:  - GA: gestational age of the fetus (in weeks)                 %
%           - motion_level: amplitude of rigid fetal movements to         %
%                           simulate                                      %
%           - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
%           - shift_mm: displacement (in mm) of the slice slab in the     %
%                       slice thickness direction between two             %
%                       simulations in the same orientation               %
%                                                                         %
%  output:  - output_folder: folder where the simulated images are saved  %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-20                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function output_folder = output_name(          GA, ...
                                     motion_level, ...
                                      orientation, ...
                                         shift_mm)

% Input check
if nargin < 4
    error('Missing input(s).');
elseif nargin > 4
    error('Too many inputs.');
end

switch motion_level
    case 0
        motion = 'no_motion';
    case 1
        motion = 'low_motion';
	case 2
        motion = 'moderate_motion';
	case 3
        motion = 'strong_motion';
end

switch orientation
    case 1
        acq_plane = 'SAG';
    case 2
        acq_plane = 'COR';
    case 3
        acq_plane = 'AX';
end

switch shift_mm
    case 0
        shift = 'no_shift';
    case -1.6
        shift = 'shift_m';
    case 1.6
        shift = 'shift_p';
end

output_folder = strcat('/data/Simu_FSE/GA_', sprintf('%02s', num2str(GA)), 'w/', motion, '/', acq_plane, '/', shift, '/');

end