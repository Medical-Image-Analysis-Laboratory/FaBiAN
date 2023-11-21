%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that implements white matter maturation processes by tuning   %
%  the T1 and T2 values from the reference maps. Two methods are          %
%  available: a 2-class segmentation using a gaussian mixture model       %
%  (GMM) and a 3-class segmentation using FAST-FSL.                       %
%                                                                         % 
%          [WM_T1map, WM_T2map] = WM_maturation(         model, ...       %
%                                                    ref_T1map, ...       %
%                                                    ref_T2map, ...       %
%                                                           GA, ...       %
%                                                       sub_id, ...       %
%                                                 segmentation, ...       %
%                                                       affine, ...       %
%                                                       SimRes, ...       %
%                                                  axcodes_reo, ...       %
%                                                        shift, ...       %
%                                                  orientation, ...       %
%                                               SamplingFactor)           %
%                                                                         %
%  inputs:  - model: anatomical model of the fetal brain                  %
%           - ref_T1map: reference T1 map of the fetal brain with the     %
%                        third dimension being the slice thickness        %
%                        direction                                        %
%           - ref_T2map: reference T2 map of the fetal brain with the     %
%                        third dimension being the slice thickness        %
%                        direction                                        %
%           - GA: gestational age of the fetus (in weeks)                 %
%           - sub_id: subject ID                                          %
%           - segmentation: 2 ways of improving white matter              %
%                           heterogeneity: 'FAST' or 'GMM'                %
%           - affine: affine orientation matrix of the anatomical model   %
%                     of the fetal brain                                  %
%           - SimRes: resolution of the 3D anatomical model of the fetal  %
%                     brain                                               %
%           - axcodes_reo: output axes codes of the reoriented 3D         %
%                          anatomical model of the fetal brain with the   %
%                          slice thickness encoded in the third           %
%                          dimension                                      %
%           - shift: displacement (in voxels) of the slice slab in the    %
%                    slice thickness direction                            %
%           - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
%           - SamplingFactor: factor by which the fetal brain volume and  %
%                             the B1 bias field have been upsampled       %
%           - InterpolationMethod: interpolation method used to assign a  %
%                                  value to every voxel of a resampled    %
%                                  3D volume. Images generated for the    %
%                                  qualitative evaluation by the          %
%                                  radiologists were simulated from       %
%                                  partial volume maps bilinearly         %
%                                  interpolated. To reduce the            %
%                                  computational burden that arises from  %
%                                  EPG simulations with many non-unique   %
%                                  combinations of [b1, T1, T2], the      %
%                                  images generated for the data          %
%                                  augmentation experiment were           %
%                                  simulated from partial volume maps     %
%                                  interpolated using a nearest-          %
%                                  -neighboor method.                     %
%                                                                         %
%  outputs: - WM_T1map: tuned T1 map of the fetal brain                   %
%           - WM_T2map: tuned T2 map of the fetal brain                   %
%                                                                         %
%                                                                         %
%  le Boeuf Andrés, 2022-03-23                                            %
%  andres.le.boeuf@estudiantat.upc.edu                                    %
%  Modified by Hélène Lajous, 2023-02-14                                  %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WM_T1map, WM_T2map] = WM_maturation(              model, ...
                                                        ref_T1map, ...
                                                        ref_T2map, ...
                                                               GA, ...
                                                           sub_id, ...
                                                     segmentation, ...
                                                           affine, ...
                                                           SimRes, ...
                                                      axcodes_reo, ...
                                                            shift, ...
                                                      orientation, ...
                                                   SamplingFactor, ...
                                              InterpolationMethod)

% Input check
if nargin < 13
    error('Missing input(s).');
elseif nargin > 13
    error('Too many inputs!');
end

tic;
% Flatten T1 and T2 3D maps 
unwrap_ref_T1map = ref_T1map(:);
unwrap_ref_T2map = ref_T2map(:);

if isequal(segmentation, 'FAST')
    %Import and set up Partial Volume Maps from FAST segmentation tool
    switch model
        case 'STA'
%             path = '.\data\atlas_fast_clustering\STA';
            path = 'C:\Users\admin\Desktop\AndrÃ©s\CODE AND DATA TO SIMULATE DIFERENT DATASETS\GMM and FAST segmentations on White matter\FAST_Clustering\STA';
            pve0 = niftiread(strcat(path, sprintf('%02s', num2str(GA)), '\STA', sprintf('%02s', num2str(GA)), '_WM_pve_0.nii.gz'));
            pve1 = niftiread(strcat(path, sprintf('%02s', num2str(GA)), '\STA', sprintf('%02s', num2str(GA)), '_WM_pve_1.nii.gz'));
            pve2 = niftiread(strcat(path, sprintf('%02s', num2str(GA)), '\STA', sprintf('%02s', num2str(GA)), '_WM_pve_2.nii.gz'));
        case 'FeTA_CHUV'
            if sub_id=='sub-709'
                path = '/data/bach/SimuHASTE/Data_processed/FeTA_challenge_2022-CHUV_testing_set_proc_resampled_1mm3/WM_segmentation_after_upsampling/';
            else
                path = '/data/bach/SimuHASTE/Data_processed/FeTA_challenge_2022-CHUV_testing_set_proc/WM_segmentation_after_upsampling/';
            end
            pve0 = niftiread(strcat(path, sub_id, '/', sub_id, '_WM_pve_0.nii.gz'));
            pve1 = niftiread(strcat(path, sub_id, '/', sub_id, '_WM_pve_1.nii.gz'));
            pve2 = niftiread(strcat(path, sub_id, '/', sub_id, '_WM_pve_2.nii.gz'));
        case 'FeTA'
            path = '/data/bach/SimuHASTE/Data_processed/FeTA2021_Release1and2Corrected_v4_proc/';
            pve0 = niftiread(strcat(path, sub_id, '/', sub_id, '_WM_pve_0.nii.gz'));
            pve1 = niftiread(strcat(path, sub_id, '/', sub_id, '_WM_pve_1.nii.gz'));
            pve2 = niftiread(strcat(path, sub_id, '/', sub_id, '_WM_pve_2.nii.gz'));
    end

    % Reorient the partial volume maps extracted from the WM mask of the
    % fetal brain so that the slice thickness direction is encoded in the
    % 3. dimension
    [pve0_reo, ~] = reorient_volume(   pve0, ...
                                     affine, ...
                                     SimRes, ...
                                axcodes_reo);
    [pve1_reo, ~] = reorient_volume(   pve1, ...
                                     affine, ...
                                     SimRes, ...
                                axcodes_reo);
    [pve2_reo, ~] = reorient_volume(   pve2, ...
                                     affine, ...
                                     SimRes, ...
                                axcodes_reo);
    
    % As for the 3D anatomical model of the fetal brain in the main.m,
    % shift the FOV of the partial volume maps extracted from the WM mask
    % in the slice thickness direction
    pve0_reo = FOV_shift(pve0_reo, shift, orientation);
    pve1_reo = FOV_shift(pve1_reo, shift, orientation);
    pve2_reo = FOV_shift(pve2_reo, shift, orientation);

    % As for the 3D anatomical model of the fetal brain in the main.m,
    % shift the FOV of the partial volume maps extracted from the WM mask
    % in the slice thickness direction
    pve0_reo = FOV_shift(pve0_reo, shift, orientation);
    pve1_reo = FOV_shift(pve1_reo, shift, orientation);
    pve2_reo = FOV_shift(pve2_reo, shift, orientation);
   
    % As for the 3D anatomical model of the fetal brain in the main.m,
    % up-sample the partial volume maps extracted from the WM mask in the
    % slice thickness direction
    pve0_reo = sampling_OoP(           pve0_reo, ...
                                 SamplingFactor, ...
                            InterpolationMethod);
    pve1_reo = sampling_OoP(           pve1_reo, ...
                                 SamplingFactor, ...
                            InterpolationMethod);
    pve2_reo = sampling_OoP(           pve2_reo, ...
                                 SamplingFactor, ...
                            InterpolationMethod);
    
    %Build advantage/disadvantage values
    pos_reward = pve0_reo(:) - pve1_reo(:);
    neg_reward = pve1_reo(:) - pve2_reo(:);
    
    pos_reward(pos_reward<0) = 0;
    neg_reward(neg_reward>0) = 0;
    
    %Out boundaries control of the advantage/disadvantage values
    ClipValue = 'adapt';
    pos_reward = clip_function(pos_reward, ClipValue, GA);
    neg_reward = clip_function(neg_reward, ClipValue, GA);
    
    fprintf('Applying T1 and T2 values variation through the white matter...\n');
    %cherche valeurs diffÃ©rentes de 1 pour appliquer variations : mettre
    %l'Ã©tape de +1 de la clip_function ici- On ajoute 1 pcq multiplication
    %aprÃ¨s
    pos_ind = find(pos_reward~=1);
    neg_ind = find(neg_reward~=1);
    
    % Wheight T1 and T2 reference values with the new advantage/disadvantage
    % values
    % for computational purpose, don't consider voxels where no difference
    % with the mean value
    for i=1:length(pos_ind)
        unwrap_ref_T1map(pos_ind(i)) = unwrap_ref_T1map(pos_ind(i)) * pos_reward(pos_ind(i));
        unwrap_ref_T2map(pos_ind(i)) = unwrap_ref_T2map(pos_ind(i)) * pos_reward(pos_ind(i));
    end
    
    for i=1:length(neg_ind)
        unwrap_ref_T1map(neg_ind(i)) = unwrap_ref_T1map(neg_ind(i)) * neg_reward(neg_ind(i));
        unwrap_ref_T2map(neg_ind(i)) = unwrap_ref_T2map(neg_ind(i)) * neg_reward(neg_ind(i));
    end
 elseif isequal(segmentation, 'GMM')
     %Import and set up Partial Volume Maps from GMM segmentation 
     path = '.\data\atlas_gmm_clustering';
     WM_hyd = niftiread(strcat(path, '\STA', sprintf('%02s', num2str(GA)), '\STA', sprintf('%02s', num2str(GA)), '_map1.nii'));
     WM_non_hyd = niftiread(strcat(path, '\STA', sprintf('%02s', num2str(GA)), '\STA', sprintf('%02s', num2str(GA)), '_map2.nii'));

     % We convert the unsigned integer values to doubles (uint16 --> maxvalue =
     % 2`16 - 1 = 65535
     WM_hyd = double(WM_hyd(:)) ./ (2^16 - 1);
     WM_non_hyd = double(WM_non_hyd(:)) ./ (2^16 - 1);

     %Searching for non background voxels
     pos_non_background = find((WM_hyd+WM_non_hyd)~=0);

     %if >0: T1 and T2 values will increase, otherwise decrease (as the WM is more dense, therefore "darker" white matter)
     reward_map = zeros(size(pos_non_background));
     for i=1:length(reward_map)
         reward_map(i) = WM_hyd(pos_non_background(i)) - WM_non_hyd(pos_non_background(i));
     end

     %We clip the change values to don't get outboundaries T1/T2 values
     ClipValue = 'adapt';
     reward_map = clip_function(reward_map, ClipValue, 'sigmoid', GA);

     fprintf('Applying T1 and T2 values variation through the white matter...\n');
    
     for i=1:length(reward_map)
         unwrap_ref_T1map(pos_non_background(i)) = unwrap_ref_T1map(pos_non_background(i)) * reward_map(i);
         unwrap_ref_T2map(pos_non_background(i)) = unwrap_ref_T2map(pos_non_background(i)) * reward_map(i);
     end
end

WM_T1map = reshape(unwrap_ref_T1map, size(ref_T1map));
WM_T2map = reshape(unwrap_ref_T2map, size(ref_T2map));
fprintf('Temporal Resolution improved in %0.5f seconds ! :D\n', toc);

end