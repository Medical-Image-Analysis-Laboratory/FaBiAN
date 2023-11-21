%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that converts segmented anatomical MR images of the fetal     %
%  brain to the corresponding reference T1 and T2 maps.                   %
%                                                                         %
%      [ref_T1map, ref_T2map] = tissue_to_MR(        fetal_model, ...     %
%                                                    Fetal_Brain, ...     %
%                                            Fetal_Brain_Tissues, ...     %
%                                                             GA, ...     %
%                                                         sub_id, ...     %
%                                               WM_heterogeneity, ...     %
%                                                         affine, ...     %
%                                                         SimRes, ...     %
%                                                    axcodes_reo, ...     %
%                                                          shift, ...     %
%                                                    orientation, ...     %
%                                                 SamplingFactor, ...     %
%                                                          T1_WM, ...     %
%                                                          T2_WM, ...     %
%                                                          T1_GM, ...     %
%                                                          T2_GM, ...     %
%                                                         T1_CSF, ...     %
%                                                         T2_CSF)         %
%                                                                         %
%  inputs:  - fetal_model: anatomical model of the fetal brain            %
%           - Fetal_Brain: segmented high-resolution 3D volume of the     %
%                          fetal brain                                    %
%           - Fetal_Brain_Tissues: list of tissues in the fetal brain     %
%                                  that were segmented and labeled        %
%           - GA: gestational age of the fetus (in weeks)                 %
%           - sub_id: subject ID                                          %
%           - WM_heterogeneity: boolean value that determines wether to   %
%                               implement (1) or not (0) white matter     %
%                               maturation processes in the simulated     %
%                               images                                    %
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
%           - T1_WM: T1 relaxation time of white matter at the magnetic   %
%                    field strength specified for the simulated           %
%                    acquisition                                          %
%           - T2_WM: T2 relaxation time of white matter at the magnetic   %
%                    field strength specified for the simulated           %
%                    acquisition                                          %
%           - T1_GM: T1 relaxation time of gray matter at the magnetic    %
%                    field strength specified for the simulated           %
%                    acquisition                                          %
%           - T2_GM: T2 relaxation time of gray matter at the magnetic    %
%                    field strength specified for the simulated           %
%                    acquisition                                          %
%           - T1_CSF: T1 relaxation time of cerebrospinal fluid at the    %
%                     magnetic field strength specified for the           %
%                     simulated acquisition                               %
%           - T2_CSF: T2 relaxation time of cerebrospinal fluid at the    %
%                     magnetic field strength specified for the           %
%                     simulated acquisition                               %
%                                                                         %
%  outputs: - ref_T1map: reference T1 map of the fetal brain              %
%           - ref_T2map: reference T2 map of the fetal brain              %
%                                                                         %
%  Relaxation values of the fetal brain at 1.5 T are based on the         %
%  literature:                                                            %
%  Hagmann, C. et al. T2 at MR imaging is an objective quantitative       %
%  measure of cerebral white matter signal intensity abnormality in       %
%  preterm infants at term-equivalent age. Radiology 252, 209–217 (2009)  %
%  https://doi.org/10.1148/radiol.2522080589                              %
%  Nossin-Manor, R. et al. Quantitative MRI in the very preterm brain:    %
%  Assessing tissue organization and myelination using magnetization      %
%  transfer, diffusion tensor and T1 imaging. NeuroImage 64, 505–516      %
%  (2013). https://doi.org/10.1016/j.neuroimage.2012.08.086               %
%  Blazejewska, A. et al. 3D in utero quantification of T2* relaxation    %
%  times in human fetal brain tissues for age optimized structural and    %
%  functional MRI. Magnetic Resonance in Medicine 78, 909–916 (2017).     %
%  https://doi.org/10.1002/mrm.26471                                      %
%  Yarnykh, V. et al. Quantitative assessment of normal fetal brain       %
%  myelination using fast macromolecular proton fraction mapping.         %
%  American Journal of Neuroradiology 39, 1341–1348 (2018).               %
%  https://doi.org/10.3174/ajnr.A5668                                     %
%  Vasylechko, S. et al. T2* relaxometry of fetal brain at 1.5 Tesla      %
%  using a motion tolerant method. Magnetic Resonance in Medicine 73,     %
%  1795–1802 (2015). https://doi.org/10.1002/mrm.25299                    %
%                                                                         %
%  Relaxation values of the fetal brain at 3 T are estimated from the     %
%  following references:                                                  %
%  Stanisz, G. J. et al. T1, T2 relaxation and magnetization transfer in  %
%  tissue at 3T. Magnetic Resonance in Medicine 54, 507–512 (2005).       %
%  https://doi.org/10.1002/mrm.20605                                      %
%  Rooney, W. D. et al. Magnetic field and tissue dependencies of human   %
%  brain longitudinal 1H2O relaxation in vivo. Magnetic Resonance in      %
%  Medicine 57, 308–318 (2007). https://doi.org/10.1002/mrm.21122         %
%  Shin, W. et al. Fast high-resolution T1 mapping using inversion-       %
%  -recovery look-locker echo-planar imaging at steady state:             %
%  Optimization for accuracy and reliability. Magnetic Resonance in       %
%  Medicine 61, 899–906 (2009). https://doi.org/10.1002/mrm.21836         %
%  Bojorquez, J. Z. et al. What are normal relaxation times of tissues    %
%  at 3 T? Magnetic Resonance Imaging 35, 69–80 (2017).                   %
%  http://dx.doi.org/10.1016/j.mri.2016.08.021                            %
%  Daoust, A. et al. Transverse relaxation of cerebrospinal fluid         %
%  depends on glucose concentration. Magnetic Resonance Imaging 44,       %
%  72–81 (2017). https://doi.org/10.1016/j.mri.2017.08.001                %
%  http://mriquestions.com/bo-effect-on-t1--t2.html                       %
%                                                                         %
%  The classification of fetal brain structures as grey matter (GM),      %
%  white matter (WM), or cerebrospinal fluid (CSF), was confirmed by two  %
%  experienced neuroradiologists at CHUV, Dre. Mériam Koob, who is also   %
%  expert in pediatric radiology, and Dr. Vincent Dunet.                  %
%      37 Hippocampus_L             => GM                                 %
%      38 Hippocampus_R             => GM                                 %
%      41 Amygdala_L                => GM                                 %
%      42 Amygdala_R                => GM                                 %
%      71 Caudate_L                 => GM                                 %
%      72 Caudate_R                 => GM                                 %
%      73 Putamen_L                 => GM                                 %
%      74 Putamen_R                 => GM                                 %
%      77 Thalamus_L                => GM                                 %
%      78 Thalamus_R                => GM                                 %
%      91 CorpusCallosum            => WM                                 %
%      92 Lateral_Ventricle_L       => CSF                                %
%      93 Lateral_Ventricle_R       => CSF                                %
%      94 Midbrain_L                => WM                                 %
%     100 Cerebellum_L              => WM                                 %
%     101 Cerebellum_R              => WM                                 %
%     108 Subthalamic_Nuc_L         => GM                                 %
%     109 Subthalamic_Nuc_R         => GM                                 %
%     110 Hippocampal_Comm          => WM                                 %
%     111 Fornix                    => WM                                 %
%     112 Cortical_Plate_L          => GM                                 %
%     113 Cortical_Plate_R          => GM                                 %
%     114 Subplate_L                => WM                                 %
%     115 Subplate_R                => WM                                 %
%     116 Inter_Zone_L              => WM                                 %
%     117 Inter_Zone_R              => WM                                 %
%     118 Vent_Zone_L (ependyme)    => WM                                 %
%     119 Vent_Zone_R               => WM                                 %
%     120 White_Matter_L            => WM                                 %
%     121 White_Matter_R            => WM                                 %
%     122 Internal_Capsule_L        => WM                                 %
%     123 Internal_Capsule_R        => WM                                 %
%     124 CSF                       => CSF                                %
%     125 misc                      => WM (as most of this label          %
%                                          corresponds to the internal    %
%                                          capsule)                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2020-05-10                                              %
%  Adapted from: XCAT_to_MR.m (https://github.com/cwroy/Fetal-XCMR/)      %
%  helene.lajous@unil.ch                                                  %
%  Modified by Andrés le Boeuf, 2022-03-23                                %
%  andres.le.boeuf@estudiantat.ucp.edu                                    %
%  Modified by Hélène Lajous, 2023-03-23                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ref_T1map, ref_T2map] = tissue_to_MR(        fetal_model, ...
                                                       Fetal_Brain, ...
                                               Fetal_Brain_Tissues, ...
                                                                GA, ...
                                                            sub_id, ...
                                                  WM_heterogeneity, ...
                                                            affine, ...
                                                            SimRes, ...
                                                       axcodes_reo, ...
                                                             shift, ...
                                                       orientation, ...
                                                    SamplingFactor, ...
                                               InterpolationMethod, ...
                                                             T1_WM, ...
                                                             T2_WM, ...
                                                             T1_GM, ...
                                                             T2_GM, ...
                                                            T1_CSF, ...
                                                            T2_CSF)

% Input check
if nargin < 19
    error('Missing input(s).');
elseif nargin > 19
    error('Too many inputs!');
end

% Assign T1 and T2 values to every fetal tissue according to its label
switch fetal_model
    case 'STA'
        Tissue(37,:) = [T1_GM, T2_GM];
        Tissue(38,:) = [T1_GM, T2_GM];
        Tissue(41,:) = [T1_GM, T2_GM];
        Tissue(42,:) = [T1_GM, T2_GM];
        Tissue(71,:) = [T1_GM, T2_GM];
        Tissue(72,:) = [T1_GM, T2_GM];
        Tissue(73,:) = [T1_GM, T2_GM];
        Tissue(74,:) = [T1_GM, T2_GM];
        Tissue(77,:) = [T1_GM, T2_GM];
        Tissue(78,:) = [T1_GM, T2_GM];
        Tissue(91,:) = [T1_WM, T2_WM];
        Tissue(92,:) = [T1_CSF, T2_CSF];
        Tissue(93,:) = [T1_CSF, T2_CSF];
        Tissue(94,:) = [T1_WM, T2_WM];
        Tissue(100,:) = [T1_WM, T2_WM];
        Tissue(101,:) = [T1_WM, T2_WM];
        Tissue(108,:) = [T1_GM, T2_GM];
        Tissue(109,:) = [T1_GM, T2_GM];
        Tissue(110,:) = [T1_WM, T2_WM];
        Tissue(111,:) = [T1_WM, T2_WM];
        Tissue(112,:) = [T1_GM, T2_GM];
        Tissue(113,:) = [T1_GM, T2_GM];
        Tissue(114,:) = [T1_WM, T2_WM];
        Tissue(115,:) = [T1_WM, T2_WM];
        Tissue(116,:) = [T1_WM, T2_WM];
        Tissue(117,:) = [T1_WM, T2_WM];
        Tissue(118,:) = [T1_WM, T2_WM];
        Tissue(119,:) = [T1_WM, T2_WM];
        Tissue(120,:) = [T1_WM, T2_WM];
        Tissue(121,:) = [T1_WM, T2_WM];
        Tissue(122,:) = [T1_WM, T2_WM];
        Tissue(123,:) = [T1_WM, T2_WM];
        Tissue(124,:) = [T1_CSF, T2_CSF];
        Tissue(125,:) = [T1_WM, T2_WM];
    case 'FeTA_CHUV'
        Tissue(1,:) = [T1_CSF, T2_CSF]; %CSF
        Tissue(2,:) = [T1_GM, T2_GM];   %GM
        Tissue(3,:) = [T1_WM, T2_WM];   %WM
        Tissue(4,:) = [T1_CSF, T2_CSF]; %lateral ventricles
        Tissue(5,:) = [T1_WM, T2_WM];   %cerebellum
        Tissue(6,:) = [T1_GM, T2_GM];   %subcortical GM
        Tissue(7,:) = [T1_WM, T2_WM];   %brainstem
    case 'FeTA' %refined FeTA dataset (Lucas Fidon, FeTA2021_Release1and2Corrected_v4)
        Tissue(1,:) = [T1_WM, T2_WM];   %WM (excluding corpus callosum)
        Tissue(2,:) = [T1_CSF, T2_CSF]; %intra-axial CSF
        Tissue(3,:) = [T1_WM, T2_WM];   %cerebellum
        Tissue(4,:) = [T1_CSF, T2_CSF]; %extra-axial CSF
        Tissue(5,:) = [T1_GM, T2_GM];   %cortical GM
        Tissue(6,:) = [T1_GM, T2_GM];   %deep GM
        Tissue(7,:) = [T1_WM, T2_WM];   %brainstem
        Tissue(8,:) = [T1_WM, T2_WM];   %corpus callosum
end

% Computation time
tic

% Initialize reference T1 and T2 maps
Unwrap_Fetal_Brain = Fetal_Brain(:);
ref_T1map = zeros(length(Unwrap_Fetal_Brain),1);
ref_T2map = zeros(length(Unwrap_Fetal_Brain),1);

% % Relaxometry properties depend on the magnetic field strength
% switch B0
%     case 1.5
%         T2_index = 2;
%         T1_index = 1;
%     case 3
%         T2_index = 4;
%         T1_index = 3;
% end

% Fetal brain properties
for label=Fetal_Brain_Tissues(1:length(Fetal_Brain_Tissues))
    brain = find(Unwrap_Fetal_Brain==label);
    for k=1:length(brain)
        ref_T1map(brain(k)) = Tissue(label, 1);
        ref_T2map(brain(k)) = Tissue(label, 2);
    end
end

ref_T1map = reshape(ref_T1map, size(Fetal_Brain));
ref_T2map = reshape(ref_T2map, size(Fetal_Brain));

% White matter maturation processes implementation
if WM_heterogeneity == 1
    [ref_T1map, ref_T2map] = WM_maturation(        fetal_model, ...
                                                     ref_T1map, ...
                                                     ref_T2map, ...
                                                            GA, ...
                                                        sub_id, ...
                                                        'FAST', ...
                                                        affine, ...
                                                        SimRes, ...
                                                   axcodes_reo, ...
                                                         shift, ...
                                                   orientation, ...
                                                SamplingFactor, ...
                                           InterpolationMethod);
end

% Display computation time
time1 = toc;
fprintf('Computation time to convert segmented high-resolution anatomical images of the fetal brain to MR contrast: %0.5f seconds.\n', time1);

end