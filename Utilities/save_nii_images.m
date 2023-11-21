%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that saves the T2w simulation of the fetal brain along with   %
%  its label map according to the gestational age of the fetus and to     %
%  the acquisition plane                                                  %
%                                                                         %
%                save_nii_images(FSE_Images, ...                          %
%                           Label_Map_recon, ...                          %
%                                        GA, ...                          %
%                               orientation, ...                          %
%                             output_folder)                              %
%                                                                         %
%  inputs:  - FSE_images: simulated T2-weighted MR images of the fetal    %
%                         brain based on the acquisition scheme of FSE    %
%                         sequences                                       %
%           - Label_Map_recon: Simulation's label map derived from the    %
%                              reference model                            %
%           - GA: gestational age of the fetus (in weeks)                 %
%           - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
%           - output_folder: folder where the simulated images are saved  %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2022-12-02                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_nii_images(       FSE_Images, ...
                           Label_Map_recon, ...
                          Simu_Binary_Mask, ...
                               orientation, ...
                                   FOVRead, ...
                                  FOVPhase, ...
                         PhaseOversampling, ...
                               reconMatrix, ...
                            SliceThickness, ...
                                  SliceGap, ...
                             output_folder, ...
                                    sub_id, ...
                                    ses_id, ...
                                    run_id)

% Input check
if nargin < 14
    error('Missing input(s).');
elseif nargin > 14
    error('Too many inputs.');
end

switch orientation
    case 1
        acq_plane = 'sagittal';
    case 2
        acq_plane = 'coronal';
    case 3
        acq_plane = 'axial';
end

% Histogram normalization between 0 and 255 of simulated SS-FSE images
norm_HASTE_Images = FSE_Images * 255 / max(max(max(abs(FSE_Images))));

% Only keep the magnitude image
HASTE_im = abs(norm_HASTE_Images);

output_path = strcat(output_folder, char(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/');
if not(isfolder(output_path))
    mkdir(output_path)
end
derivatives_path = strcat(output_folder, 'derivatives/');
if not(isfolder(derivatives_path))
    mkdir(derivatives_path)
end
derivatives_labels_path = strcat(derivatives_path, 'labels/', char(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/');
derivatives_masks_path = strcat(derivatives_path, 'masks/', char(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/');
if not(isfolder(derivatives_labels_path))
    mkdir(derivatives_labels_path)
end
if not(isfolder(derivatives_masks_path))
    mkdir(derivatives_masks_path)
end

output_im = strcat(output_path, char(sub_id), '_ses-', sprintf('%02s', num2str(ses_id)), '_run-', num2str(run_id), '_T2w.nii');
output_labels = strcat(derivatives_labels_path, char(sub_id), '_ses-', sprintf('%02s', num2str(ses_id)), '_run-', num2str(run_id), '_labels.nii');
output_mask = strcat(derivatives_masks_path, char(sub_id), '_ses-', sprintf('%02s', num2str(ses_id)), '_run-', num2str(run_id), '_mask.nii');

% Save images
gzip(mat2nii(HASTE_im, output_im, acq_plane, [FOVRead/reconMatrix (FOVPhase/(1+PhaseOversampling))/reconMatrix SliceThickness+SliceGap]));
gzip(mat2nii(Label_Map_recon, output_labels, acq_plane, [FOVRead/reconMatrix (FOVPhase/(1+PhaseOversampling))/reconMatrix SliceThickness+SliceGap]));
gzip(mat2nii(Simu_Binary_Mask, output_mask, acq_plane, [FOVRead/reconMatrix (FOVPhase/(1+PhaseOversampling))/reconMatrix SliceThickness+SliceGap]));

end