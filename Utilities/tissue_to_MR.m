%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that converts segmented anatomical MR images of the fetal     %
%  brain to the corresponding reference T1 and T2 maps depending on the   %
%  main magnetic field strength B0.                                       %
%                                                                         %
%       [ref_T1map, ref_T2map] = tissue_to_MR(        Fetal_Brain, ...    %
%                                             Fetal_Brain_Tissues, ...    %
%                                                              B0);       %
%                                                                         %
%  input:   - Fetal_Brain: segmented high-resolution 3D volume of the     %
%                          fetal brain                                    %
%           - Fetal_Brain_Tissues: list of tissues in the fetal brain     %
%                                  that were segmented and labeled        %
%           - B0: main magnetic field strength                            %
%                                                                         %
%  outputs: - ref_T1map: reference T1 map of the fetal brain              %
%           - ref_T2map: reference T2 map of the fetal brain              %
%                                                                         %
%  T1 and T2 values of the segmented fetal brain tissues are stored in a  %
%  table with the format:                                                 %
%  [T1 at 1.5 T, T2 at 1.5 T, T1 at 3.0 T, T2 at 3.0 T]                   %
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
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ref_T2map, ref_T1map] = tissue_to_MR(        Fetal_Brain, ...
                                               Fetal_Brain_Tissues, ...
                                                                B0)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs!');
end

% Assign T1 and T2 values to every fetal tissue according to its label
Tissue(37,:) = [2000, 162, 2500, 162];
Tissue(38,:) = [2000, 162, 2500, 162];
Tissue(41,:) = [2000, 162, 2500, 162];
Tissue(42,:) = [2000, 162, 2500, 162];
Tissue(71,:) = [2000, 162, 2500, 162];
Tissue(72,:) = [2000, 162, 2500, 162];
Tissue(73,:) = [2000, 162, 2500, 162];
Tissue(74,:) = [2000, 162, 2500, 162];
Tissue(77,:) = [2000, 162, 2500, 162];
Tissue(78,:) = [2000, 162, 2500, 162];
Tissue(91,:) = [3000, 232, 3300, 232];
Tissue(92,:) = [4000, 2000, 4400, 2000];
Tissue(93,:) = [4000, 2000, 4400, 2000];
Tissue(94,:) = [3000, 232, 3300, 232];
Tissue(100,:) = [3000, 232, 3300, 232];
Tissue(101,:) = [3000, 232, 3300, 232];
Tissue(108,:) = [2000, 162, 2500, 162];
Tissue(109,:) = [2000, 162, 2500, 162];
Tissue(110,:) = [3000, 232, 3300, 232];
Tissue(111,:) = [3000, 232, 3300, 232];
Tissue(112,:) = [2000, 162, 2500, 162];
Tissue(113,:) = [2000, 162, 2500, 162];
Tissue(114,:) = [3000, 232, 3300, 232];
Tissue(115,:) = [3000, 232, 3300, 232];
Tissue(116,:) = [3000, 232, 3300, 232];
Tissue(117,:) = [3000, 232, 3300, 232];
Tissue(118,:) = [3000, 232, 3300, 232];
Tissue(119,:) = [3000, 232, 3300, 232];
Tissue(120,:) = [3000, 232, 3300, 232];
Tissue(121,:) = [3000, 232, 3300, 232];
Tissue(122,:) = [3000, 232, 3300, 232];
Tissue(123,:) = [3000, 232, 3300, 232];
Tissue(124,:) = [4000, 2000, 4400, 2000];
Tissue(125,:) = [3000, 232, 3300, 232];

% Computation time
tic

% Initialize reference T1 and T2 maps
Unwrap_Fetal_Brain = Fetal_Brain(:);
ref_T1map = zeros(length(Unwrap_Fetal_Brain),1);
ref_T2map = zeros(length(Unwrap_Fetal_Brain),1);

% Relaxometry properties depend on the magnetic field strength
switch B0
    case 1.5
        T2_index = 2;
        T1_index = 1;
    case 3
        T2_index = 4;
        T1_index = 3;
end

% Fetal brain properties
for label=Fetal_Brain_Tissues(1:length(Fetal_Brain_Tissues))
    brain = find(Unwrap_Fetal_Brain==label);
    for k=1:length(brain)
        ref_T1map(brain(k)) = Tissue(label, T1_index);
        ref_T2map(brain(k)) = Tissue(label, T2_index);
    end
end
ref_T1map = reshape(ref_T1map, size(Fetal_Brain));
ref_T2map = reshape(ref_T2map, size(Fetal_Brain));

% Display computation time
fprintf('Computation time to convert segmented high-resolution anatomical images of the fetal brain to MR contrast: %0.5f seconds.\n', toc);

end