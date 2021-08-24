%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that returns all information on the rigid motion experienced  %
%  by a 3D object depending on the amplitude of motion specified by the   %
%  user.                                                                  %
%                                                                         %
%       [ motion_corrupted_slices, ...                                    %
%        translation_displacement, ...                                    %
%                  rotation_angle, ...                                    %
%                   rotation_axis] = motion_transform(motion_level, ...   %
%                                                         NbSlices)       %
%                                                                         %
%  input:   - motion_level: amplitude of rigid motion of the fetus to be  %
%                           simulated (0: no motion; 1: little motion;    %
%                           2: moderate motion; 3: strong motion)         %
%           - NbSlices: number of slices                                  %
%                                                                         %
%  output:  - motion_corrupted_slices: list of indexes of the slices      %
%                                      corrupted by motion                %
%           - translation_displacement: table of translation              %
%                                       displacements along the 3 main    %
%                                       axes for each motion-corrupted    %
%                                       slice                             %
%           - rotation_angle: list of rotation angles for each            %
%                             motion-corrupted slice                      %
%           - rotation_axis: table of rotation axes for each              %
%                            motion-corrupted slice                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ motion_corrupted_slices, ...
          translation_displacement, ...
                    rotation_angle, ...
                     rotation_axis] = motion_transform(motion_level, ...
                                                           NbSlices)


% Input check
if nargin < 2
    error('Missing input(s).');
elseif nargin > 2
    error('Too many inputs.');
end


% We define 3 levels of motion:
%  1 - Slight motion of the fetus: in this case, less than 5% of slices are
%      affected by fetal motion, translation will vary between [-1,+1]mm
%      and rotation between [-2,+2]°;
%  2 - Moderate motion: less than 10% of slices are affected, translation
%  will vary between [-2,+2]mm and rotation between [-4,+4]°;
%  3 - Strong motion: less than 25% of slices are affected, translation
%  will vary between [-4,+4]mm and rotation between [-8,+8]°.
switch motion_level
    case 0  %no motion
        slice_percentage = 0;
    case 1  %slight motion
        slice_percentage = ceil(NbSlices*5/100);
        translation_amplitude = 1;
        rotation_amplitude = 2;
    case 2  %moderate motion
        slice_percentage = ceil(NbSlices*5/100);
        translation_amplitude = 3;
        rotation_amplitude = 5;
    case 3  %strong motion
        slice_percentage = ceil(NbSlices*5/100);
        translation_amplitude = 4;
        rotation_amplitude = 8;
end

% List of slices corrupted by motion
motion_corrupted_slices = randperm(NbSlices, slice_percentage(randperm(length(slice_percentage),1)));

% Initialization of tranformation arrays
translation_displacement = zeros(length(motion_corrupted_slices), 3);
rotation_axis = zeros(length(motion_corrupted_slices), 3);
rotation_angle = zeros(length(motion_corrupted_slices), 1);

for iSlice=1:length(motion_corrupted_slices)
    % List of translation displacement along the 3 main axes for each
    % motion-corrupted slice
    translation_displacement(iSlice,:) = [unifrnd(-translation_amplitude,translation_amplitude) ...
                                          unifrnd(-translation_amplitude,translation_amplitude) ...
                                          unifrnd(-translation_amplitude,translation_amplitude)];    %in mm
    % Rotation
    rotation_angle(iSlice) = unifrnd(-rotation_amplitude,rotation_amplitude);   %in degrees
    rotation_axis(iSlice,:) = [unifrnd(0,1) unifrnd(0,1) unifrnd(0,1)];
    % Display message for debugging
    sprintf('The motion corrupted slice #: %d is characterized by a translational displacement of [%0.4f %0.4f %0.4f]mm along the 3 main axes, and by a rotation of %0.4f degrees around the axis [%0.4f %0.4f %0.4f].', motion_corrupted_slices(iSlice), translation_displacement(iSlice,:), rotation_angle(iSlice), rotation_axis(iSlice,:))
end
