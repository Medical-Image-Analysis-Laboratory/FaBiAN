%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This demonstration script of FaBiAN use and functionalities generates  %
%  T2-weighted MR images of the fetal brain based on highly flexible      %
%  simulations of fast spin echo (FSE) sequences from various MR vendors  %
%  and at different magnetic field strengths throughout gestation.        %
%  This demonstration script simulates:                                   %
%  (i) HASTE acquisitions as implemented by Siemens Healthineers on a     %
%  1.5-T MAGNETOM Sola system (Siemens Healthineers, Erlangen, Germany)   %
%  at Lausanne University Hospital (CHUV);                                %
%  (ii) SS-FSE sequences as implemented by GE Healthcare on a 3-T whole-  %
%  -body scanner (Signa Discovery MR750) at University Children’s         %
%  Hospital Zurich.                                                       %
%                                                                         %
%  Images are simulated for a fetus of gestational age (GA):              %
%  (i) 26 weeks in the axial orientation without any shift of the slice   %
%  slab, and with random little motion;                                   %
%  (ii) 33 weeks in the sagittal orientation with a shift of -1.6 mm of   %
%  the slice slab, and with strong motion.                                %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-07-19                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  (i) HASTE acquisition (Siemens Healthineers)                           %
%      B0=1.5T;                                                           %
%      GA=26weeks;                                                        %
%      orientation=axial; no shift                                        %
%      little motion of the fetus                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

addpath('Utilities')

% Fetal brain model: In this demonstration, we base our simulations on
% segmented high-resolution anatomical MR images of the fetal brain that
% can be downloaded from:
% http://crl.med.harvard.edu/research/fetal_brain_atlas/
% Gholipour, A. et al. A normative spatiotemporal MRI atlas of the fetal
% brain for automatic segmentation and analysis of early brain growth.
% Scientific Reports 7, 476 (2017).
% https://doi.org/10.1038/s41598-017-00525-w
Fetal_Brain_model_path = 'data/Simu_FSE/Atlas/CRL_FetalBrainAtlas_2017v3/';
% Gestational age (in weeks)
GA = 26;
% Resolution of the Fetal_Brain images (isotropic, in mm)
SimRes = 0.8;
% Introduce a shift variable to slightly shift the slice series between two
% simulations in the same orientation
shift_mm = 0;   %mm
% Choose the orientation plane of the acquisitions
% (1: sagittal, 2: coronal, 3: axial)
orientation = 3;
% Non-linear slowly-varying intensity non-uniformity (INU) fields (b1+) can
% be downloaded from BrainWeb database:
% https://brainweb.bic.mni.mcgill.ca/brainweb/about_sbd.html
inu = 'data/Simu_FSE/rf20_B.rawb';
% Define a sampling factor to subdivide the volume in the slice thickness
% orientation
sampling_factor = 1;
% Main magnetic field strength
B0 = 1.5;
% Acquisition parameters
ESP = 4.08;  %ms
ETL = 224;
% Geometry
PhaseOversampling = 0.803571000000000;
SliceThickness = 3.2; %mm
SliceGap = 0.3; %mm
% Resolution
FOVRead = 360;  %mm
FOVPhase = 360; %mm
BaseResolution = 327;   %voxels
PhaseResolution = 0.7;
% Contrast
TR = 4.08;  %ms
TEeff = 90; %ms
% Acceleration technique
ACF = 2;
RefLines = 42;
% Motion
motion_level = 3;   %little motion
% Scanner zero-interpolation filling (ZIP)
% (0: no ZIP; 1: Fermi filtering in k-space and ZIP)
zip = 0;
reconMatrix = BaseResolution;
% SNR
std_noise = 0;

output_folder = output_name(          GA, ...
                            motion_level, ...
                             orientation, ...
                                shift_mm);
cd ..

mkdir(output_folder)
cd FaBiAN
t1=cputime;
HASTE_Images = FaBiAN_main(Fetal_Brain_model_path, ...
                                               GA, ...
                                           SimRes, ...
                                         shift_mm, ...
                                      orientation, ...
                                              inu, ...
                                  sampling_factor, ...
                                               B0, ...
                                              ESP, ...
                                              ETL, ...
                                PhaseOversampling, ...
                                   SliceThickness, ...
                                         SliceGap, ...
                                          FOVRead, ...
                                         FOVPhase, ...
                                   BaseResolution, ...
                                  PhaseResolution, ...
                                               TR, ...
                                            TEeff, ...
                                              ACF, ...
                                         RefLines, ...
                                     motion_level, ...
                                              zip, ...
                                      reconMatrix, ...
                                        std_noise, ...
                                    output_folder);
tfinal=cputime-t1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  (ii) SS-FSE sequence (GE Healthcare)                                   %
%       B0=3T;                                                            %
%       GA=33weeks;                                                       %
%       orientation=sagittal; shift=-1.6mm                                %
%       strong motion of the fetus                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

addpath('Utilities')

% Fetal brain model: In this demonstration, we base our simulations on
% segmented high-resolution anatomical MR images of the fetal brain that
% can be downloaded from:
% http://crl.med.harvard.edu/research/fetal_brain_atlas/
% Gholipour, A. et al. A normative spatiotemporal MRI atlas of the fetal
% brain for automatic segmentation and analysis of early brain growth.
% Scientific Reports 7, 476 (2017).
% https://doi.org/10.1038/s41598-017-00525-w
Fetal_Brain_model_path = 'data/Simu_FSE/Atlas/CRL_FetalBrainAtlas_2017v3/';
% Gestational age (in weeks)
GA = 33;
% Resolution of the Fetal_Brain images (isotropic, in mm)
SimRes = 0.8;
% Introduce a shift variable to slightly shift the FOV between 2
% simulations
shift_mm = -1.6;    %mm
% Choose the orientation plane of the acquisitions
% (1: sagittal, 2: coronal, 3: axial)
orientation = 1;
% Non-linear slowly-varying intensity non-uniformity (INU) fields (b1+) can
% be downloaded from BrainWeb database:
% https://brainweb.bic.mni.mcgill.ca/brainweb/about_sbd.html
inu = 'data/Simu_FSE/rf20_B.rawb';
% Define a sampling factor to subdivide the volume in the slice thickness
% orientation
sampling_factor = SimRes/0.1;
% Main magnetic field strength
B0 = 3;
% Acquisition parameters
ESP = 10;  %ms
ETL = 224;
% Geometry
PhaseOversampling = 0;
SliceThickness = 3; %mm
SliceGap = 0; %mm
% Resolution
FOVRead = 240.0256;  %mm
FOVPhase = 240.0256; %mm
BaseResolution = 256;   %voxels
PhaseResolution = 1;
% Contrast
TR = 10;  %ms
TEeff = 118.08; %ms
% Acceleration technique
ACF = 1;
RefLines = 0;
% Motion
motion_level = 3;
% Scanner zero-interpolation filling (ZIP)
% (0: no ZIP; 1: Fermi filtering in k-space and ZIP)
zip = 1;
reconMatrix = BaseResolution*2;
% SNR
std_noise = 0.01;

output_folder = output_name(          GA, ...
                            motion_level, ...
                             orientation, ...
                                shift_mm);
cd FaBiAN
SSFSE_Images = FaBiAN_main(Fetal_Brain_model_path, ...
                                               GA, ...
                                           SimRes, ...
                                         shift_mm, ...
                                      orientation, ...
                                              inu, ...
                                  sampling_factor, ...
                                               B0, ...
                                              ESP, ...
                                              ETL, ...
                                PhaseOversampling, ...
                                   SliceThickness, ...
                                         SliceGap, ...
                                          FOVRead, ...
                                         FOVPhase, ...
                                   BaseResolution, ...
                                  PhaseResolution, ...
                                               TR, ...
                                            TEeff, ...
                                              ACF, ...
                                         RefLines, ...
                                     motion_level, ...
                                              zip, ...
                                      reconMatrix, ...
                                        std_noise, ...
                                    output_folder);
