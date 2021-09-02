%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that redefines the position of the slice slab that covers     %
%  the fetal brain volume based on the shift that is applied in the       %
%  slice thickness direction between the acquisition of two               %
%  low-resolution series in the same orientation (as it is done in        %
%  clinical routine in the case of two successive series acquired in the  %
%  same plane).                                                           %
%                                                                         %
%     Fetal_Brain_shift = FOV_shift(Fetal_Brain, shift, orientation);     %
%                                                                         %
%  inputs:  - Fetal_Brain: segmented high-resolution 3D volume of the     %
%                          fetal brain                                    %
%           - shift: displacement (in voxels) of the slice slab in the    %
%                    slice thickness direction                            %
%           - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
%                                                                         %
%  output:  - Fetal_Brain_shift: segmented high-resolution 3D volume of   %
%                                the fetal brain after a shift in the     %
%                                slice thickness direction                %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-08-23                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Fetal_Brain_shift = FOV_shift(Fetal_Brain, ...
                                             shift, ...
                                       orientation)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs.');
end

% Resize the fetal brain volume to allow shifting in any direction
Fetal_Brain_FOV = {1+ceil(abs(shift)):size(Fetal_Brain,1)-ceil(abs(shift)); 1+ceil(abs(shift)):size(Fetal_Brain,2)-ceil(abs(shift)); 1+ceil(abs(shift)):size(Fetal_Brain,3)-ceil(abs(shift))};

% Slightly shift the slice slab in the slice thickness direction (as done
% in the clinics)
Fetal_Brain_FOV{orientation} = Fetal_Brain_FOV{orientation} + shift;

% Mask the fetal brain after having shifted the slice slab
Fetal_Brain_shift = Fetal_Brain(Fetal_Brain_FOV{1}, Fetal_Brain_FOV{2}, Fetal_Brain_FOV{3});

% Display message for debugging
sprintf('The fetal brain volume has been shifted by %d voxels in the slice thickness direction.', shift)

end