%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that updates the orientation matrix of a modified 3D image    %
%  (e.g., after reorientation, zero-padding, etc.).                       %
%                                                                         %
%                NewAffine = update_affine(      Model, ...               % 
%                                          ModelAffine, ...               %
%                                             NewImage, ...               %
%                                           InitAffine, ...               %
%                                         CenterOffset)                   %
%                                                                         %
%  inputs:  - Model: 3D volume with know affine                           %
%           - ModelAffine: affine transformation matrix corresponding to  %
%                          the Model volume                               %
%           - NewImage: volume generated from the information of the      %
%                       Model                                             %
%           - InitAffine: initialization of the affine matrix of the      %
%                         NewImage volume. Information on the offset      %
%                         values will be computed by this function.       %
%           - CenterOffset: offset to substract to the center of the      %
%                           NewImage computed in intrinsic coordinates    %
%                           to account for any deviation due to           %
%                           zero-padding steps occuring across the        %
%                           simulation pipeline                           %
%                                                                         %
%  output:  - NewAffine: full affine transformation matrix of the         %
%                        NewImage volume.                                 %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2023-02-27                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function NewAffine = update_affine(       Model, ...
                                    ModelAffine, ...
                                       NewImage, ...
                                   NewRotations, ...
                                       varargin)

% Input check
if nargin < 4
    error('Missing input(s).');
elseif nargin==4
    CenterOffset = 0; %default value if the center of the NewImage has not been modified
else
    if numel(varargin) > 0  %optional input arguments are provided
        if numel(varargin) < 2
            error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
        end
        switch varargin{1}
            case 'CenterOffset'
                CenterOffset = varargin{2};
            otherwise
                error('Unexpected ''ParameterName'' input: %s\n', varargin{1});
        end
    end
end

% Compute the center of the original anatomical images in intrinsic
% coordinates
ModelCenter_vox = (size(Model) - 1) / 2;  %The origin of a coordinate system is (0,0,0)
ModelCenter_vox(4) = 1;

% Compute the center of the original anatomical images in world coordinates
ModelCenter_mm = ModelAffine * ModelCenter_vox';

% % Initialize the orientation matrix of the simulated images
% NewRotations = (InitAffine(1:3, 1:3) ./ SimResReo) .* SimVoxSize;

% Initialize the center of the simulated images in intrinsic coordinates
% if ndims(NewImage)>3
%     NewImage = squeeze(NewImage(:,:,:,1));
% end
SimCenter_vox = (size(NewImage) - 1) / 2 + CenterOffset;

% Compute the center of the simulated images in world coordinates
SimCenter_mm = NewRotations * SimCenter_vox';

% Offset
Offset = ModelCenter_mm(1:3) - SimCenter_mm;

% Define the orientation matrix of the simulated images to be in the same
% space as the original anatomical model
NewAffine = eye(size(ModelAffine));
NewAffine(1:3, 1:3) = NewRotations;
NewAffine(1:3, 4) = Offset;

end