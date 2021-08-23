%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that redefines the field-of-view (FOV) that covers the fetal  %
%  brain volume based on the shift that is applied in the slice           %
%  thickness direction between the acquisition of two low-resolution      %
%  series in the same orientation.                                        %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
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


% Slightly shift the FOV in the slice thickness direction (as done in the
% clinics)
Fetal_Brain_FOV = {1+ceil(abs(shift)):size(Fetal_Brain,1)-(1+ceil(abs(shift))); 1+ceil(abs(shift)):size(Fetal_Brain,2)-(1+ceil(abs(shift))); 1+ceil(abs(shift)):size(Fetal_Brain,3)-(1+ceil(abs(shift)))};
Fetal_Brain_FOV{orientation} = Fetal_Brain_FOV{orientation} + shift;

% Mask the fetal brain after having shifted the FOV
Fetal_Brain_shift = Fetal_Brain(Fetal_Brain_FOV{1}, Fetal_Brain_FOV{2}, Fetal_Brain_FOV{3});

% Display message for debugging
sprintf('The fetal brain volume has been shifted by %d voxels in the slice thickness direction.', shift)
