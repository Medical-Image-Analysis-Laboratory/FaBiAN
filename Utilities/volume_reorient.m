%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that reorients a 3D volume to have its orientation, which     %
%  corresponds to the slice thickness direction, given by the 3. index.   %
%                                                                         %
%        Volume_Reoriented = volume_reorient(Volume, orientation);        %
%                                                                         %
%  input:   - Volume: 3D volume to be reoriented                          %
%           - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
%                                                                         %
%  output:  - Volume_Reoriented: 3D volume after reorientation, i.e. its  %
%                                3. dimension corresponds to the slice    %
%                                thickness direction                      %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Volume_Reoriented = volume_reorient(     Volume, ...
                                             orientation)

% Input check
if nargin < 2
    error('Missing input(s).');
elseif nargin > 2
    error('Too many inputs.');
end

Volume_Reoriented = Volume;

% For the sake of clarity, the orientation corresponding to the slice
% thickness will always be given by the 3. index
switch orientation
    case 1   %sagittal
        Volume_Reoriented = permute(Volume, [2,3,1]);
    case 2   %coronal
        Volume_Reoriented = permute(Volume, [1,3,2]);
end

% Display message for debugging
sprintf('Volume reorientation.')

end