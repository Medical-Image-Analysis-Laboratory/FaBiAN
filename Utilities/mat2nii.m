%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that converts a matlab matrix into an image in .nifti format. %
%                                                                         %
%                 varargout = mat2nii(simulated_image, ...                %
%                                        voxSize_simu, ...                %
%                                       model_niiinfo, ...                %
%                                        out_filename)                    %
%                                                                         %
%  inputs:  - simulated_image: simulated FSE image stored as .mat file    %
%           - voxSize_simu: voxel size of the simulated image             %
%           - model_niiinfo: information (nifti header) of the original   %
%                            anatomical model from which the simulated    %
%                            image is derived                             %
%           - out_filename: output filename (including path) of the       %
%                           generated nifti image                         %
%                                                                         %
%  output: The Matlab matrix simulated_image is saved in nifti format in  %
%          out_filename.                                                  %
%                                                                         %
%                                                                         %
%  Oscar Esteban, Hélène Lajous, 2023-01-20                               %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function varargout = mat2nii(simulated_image, ...
                                voxSize_simu, ...
                               model_niiinfo, ...
                                out_filename)

% Input check
if nargin < 4
    error('Missing input(s).');
elseif nargin > 4
    error('Too many inputs.');
end

% Retrieve information on the original anatomical model from its header
orientation_matrix_model = model_niiinfo.Transform.T';

% Compute the center of the original anatomical images in intrinsic
% coordinates
center_vox_model = (model_niiinfo.ImageSize - 1) / 2;
center_vox_model(4) = 1;

% Compute the center of the original anatomical images in world coordinates
center_mm_model = orientation_matrix_model * center_vox_model';

% Initialize the orientation matrix of the simulated images
rotations_sim = sign(orientation_matrix_model(1:3, 1:3)) .* voxSize_simu;

% Compute the center of the simulated images in intrinsic coordinates
center_vox_sim = (size(simulated_image) - 1) / 2;

% Compute the center of the simulated images in world coordinates
center_mm_sim = rotations_sim * center_vox_sim';

% Offset
offset = center_mm_sim - center_mm_model(1:3);

% Define the orientation matrix of the simulated images to be in the same
% space as the original anatomical model
orientation_matrix_simu = eye(size(orientation_matrix_model));
orientation_matrix_simu(1:3, 1:3) = rotations_sim;
orientation_matrix_simu(1:3, 4) = -offset;

% Update the information of the simulated images to save them in nifti
% format with a proper nifti header

newinfo = model_niiinfo;
newinfo.Filename = out_filename;
newinfo.Transform.T = orientation_matrix_simu';
newinfo.TransformName = "Qform";
newinfo.ImageSize = size(simulated_image);
newinfo.PixelDimensions = voxSize_simu;
if isa(simulated_image, model_niiinfo.Datatype)==0
    newinfo.Datatype = class(simulated_image);
end
niftiwrite(simulated_image, out_filename, newinfo, 'Compressed', true);

% Output
varargout{1} = out_filename;

end