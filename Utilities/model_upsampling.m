%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that reads the segmented high-resolution anatomical MR        %
%  images of the fetal brain from which the simulated images will be      %
%  derived at a given gestational age (GA).                               %
%  Here, we consider segmented high-resolution images from Gholipour, A.  %
%  et al. A normative spatiotemporal MRI atlas of the fetal brain for     %
%  automatic segmentation and analysis of early brain growth. Scientific  %
%  Reports 7, 476 (2017). https://doi.org/10.1038/s41598-017-00525-w      %
%                                                                         %
%          Fetal_Brain = model_upsampling(             path, ...          %
%                                                    sub_id, ...          %
%                                         upsampling_factor)              %
%                                                                         %
%  inputs:  - path: folder where the segmented high-resolution MR images  %
%                   of the fetal brain from which our simulations are     %
%                   derived are stored                                    %
%           - sub_id: subject ID                                          %
%           - upsampling_factor: factor to upsample the anatomical model  %
%                                to the nearest decimal to reduce the     %
%                                computational burden                     %
%                                                                         %
%  output:  - Fetal_Brain: segmented high-resolution 3D volume of the     %
%                          fetal brain for subject sub_id                 %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2022-11-28                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Fetal_Brain = model_upsampling(             path, ...
                                                   sub_id, ...
                                        upsampling_factor)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs.');
end

info = niftiinfo(strcat(path, sub_id, '\anat\', sub_id, '_rec-mial_dseg.nii.gz'));
Fetal_Brain_init = niftiread(info);
Fetal_Brain = imresize3(Fetal_Brain_init, [size(Fetal_Brain_init,1)*upsampling_factor size(Fetal_Brain_init,2)*upsampling_factor size(Fetal_Brain_init,3)*upsampling_factor], 'nearest');
% info.PixelDimensions = [SimRes SimRes SimRes];
% info.ImageSize = [size(Fetal_Brain_up,1) size(Fetal_Brain_up,2) size(Fetal_Brain_up,3)];
% niftiwrite(Fetal_Brain_up, 'outbrain.nii', info);

end