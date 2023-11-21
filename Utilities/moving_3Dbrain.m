%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that translates and rotates a 3D volume.                      %
%                                                                         %
%    Fetal_Brain_rotated = moving_3Dbrain(              Fetal_Brain, ...  %
%                                                            SimRes, ...  %
%                                          translation_displacement, ...  %
%                                         translation_interpolation, ...  %
%                                                    rotation_angle, ...  %
%                                                     rotation_axis, ...  %
%                                            rotation_interpolation, ...  %
%                                                              bbox);     %
%                                                                         %
%  inputs:  - Fetal_Brain: segmented high-resolution 3D volume of the     %
%                          fetal brain                                    %
%           - SimRes: resolution of the original high-resolution fetal    %
%                     brain images (isotropic, in mm)                     %
%           - translation_displacement: table of translation              %
%                                       displacements (in mm) along the   %
%                                       three main axes for each motion-  %
%                                       -corrupted slice                  %
%           - translation_interpolation: interpolation method used to     %
%                                        assign an intensity value to     %
%                                        every voxel of the 3D volume     %
%                                        after translation                %
%           - rotation_angle: list of rotation angles (in °) for 3D       %
%                             rotation in every motion-corrupted slice    %
%           - rotation_axis: table of rotation axes for 3D rotation in    %
%                            every motion-corrupted slice                 %
%           - rotation_interpolation: interpolation method used to        %
%                                     assign an intensity value to every  %
%                                     voxel of the 3D volume after        %
%                                     rotation                            %
%           - bbox: bounding box that controls the size of the output     %
%                   volume after rotation                                 %
%                                                                         %
%  output:  - Fetal_Brain_rotated: segmented high-resolution 3D volume    %
%                                  of the fetal brain after 3D            %
%                                  translation and rotation               %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-20                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Fetal_Brain_rotated = moving_3Dbrain(              Fetal_Brain, ...
                                                                 SimRes, ...
                                               translation_displacement, ...
                                              translation_interpolation, ...
                                                         rotation_angle, ...
                                                          rotation_axis, ...
                                                 rotation_interpolation, ...
                                                                   bbox)

% Input check
if nargin < 8
    error('Missing input(s).');
elseif nargin > 8
    error('Too many inputs.');
end

% Use the information on the resolution of the anatomical model in each
% dimension to construct a spatial referencing object associated with it
Reference_Fetal_Brain_mm = imref3d(size(Fetal_Brain), SimRes(2), SimRes(1), SimRes(3));

% Move the fetal brain anatomy independently along the three main axes
% according to the given translation displacements (in mm)
[Fetal_Brain_translated, ~] = imtranslate(Fetal_Brain, Reference_Fetal_Brain_mm, translation_displacement, translation_interpolation);

% Apply 3D rotation, defined by a rotation angle and a rotation axis, to
% the already translated fetal brain anatomy
Fetal_Brain_rotated = imrotate3(Fetal_Brain_translated, rotation_angle, rotation_axis, rotation_interpolation, bbox);

end