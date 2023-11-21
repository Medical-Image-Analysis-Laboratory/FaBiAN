%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This is FaBiAN's main script to generate T2-weighted MR images of the  %
%  fetal brain by simulating the physical principles involved in fast     %
%  spin echo (FSE) sequences.                                             %
%  Our numerical framework is based on the extended phase graph (EPG)     %
%  formalism described in details in: Weigel, M. Extended phase graphs:   %
%  Dephasing, RF pulses, and echoes - pure and simple. Journal of         %
%  Magnetic Resonance Imaging 41, 266-295 (2015).                         %
%  https://doi.org/10.1002/jmri.24619. The associated EPG simulation      %
%  code from Matthias Weigel for multi-spin echo sequences can be         %
%  downloaded here: https://github.com/matthias-weigel/EPG.               %
%                                                                         %
%        function FSE_Images = FaBiAN_main(Fetal_Brain_model_path, ...    %
%                                                              GA, ...    %
%                                                          SimRes, ...    %
%                                                        shift_mm, ...    %
%                                                     orientation, ...    %
%                                                             inu, ...    %
%                                                  SamplingFactor, ...    %
%                                                              B0, ...    %
%                                                             ESP, ...    %
%                                                             ETL, ...    %
%                                               PhaseOversampling, ...    %
%                                                  SliceThickness, ...    %
%                                                        SliceGap, ...    %
%                                                         FOVRead, ...    %
%                                                        FOVPhase, ...    %
%                                                  BaseResolution, ...    %
%                                                 PhaseResolution, ...    %
%                                                              TR, ...    %
%                                                           TEeff, ...    %
%                                                             ACF, ...    %
%                                                        RefLines, ...    %
%                                                    motion_level, ...    %
%                                                             zip, ...    %
%                                                     ReconMatrix, ...    %
%                                                       std_noise, ...    %
%                                                   output_folder);       %
%                                                                         %
%  inputs:  - Fetal_Brain_model_path: folder where the segmented high-    %
%                                     -resolution MR images of the fetal  %
%                                     brain from which our simulations    %
%                                     are derived are stored              %
%           - GA: gestational age of the fetus (in weeks)                 %
%           - SimRes: resolution of the original high-resolution fetal    %
%                     brain images (isotropic, in mm)                     %
%           - shift_mm: displacement (in mm) of the slice slab in the     %
%                       slice thickness direction                         %
%           - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
%           - inu: simulated intensity non-uniformity fields (3D)         %
%           - SamplingFactor: factor of up- or down-sampling              %
%           - B0: main magnetic field strength                            %
%           - ESP: echo spacing (in ms)                                   %
%           - ETL: echo train length, i.e. number of 180Â°-RF pulses       %
%           - PhaseOversampling: oversampling in the phase-encoding       %
%                                direction                                %
%           - SliceThickness: thickness of a slice (in mm)                %
%           - SliceGap: distance (in mm) between two consecutive slices   %
%           - FOVRead: dimension of the field-of-view in the read-out     %
%                      direction (in mm)                                  %
%           - FOVPhase: dimension of the field-of-view in the phase-      %
%                       encoding direction (in mm)                        %
%           - BaseResolution: matrix size in the read-out direction (in   %
%                             voxels)                                     %
%           - PhaseResolution: matrix size in the phase-encoding          %
%                              direction (* BaseResolution)               %
%           - TR: time from the application of an excitation pulse to     %
%                 the application of the next pulse (i.e., echo spacing   %
%                 in the EPG simulations) (in ms)                         %
%           - TEeff: effective echo time (in ms)                          %
%           - ACF: acceleration factor                                    %
%           - RefLines: number of lines that are consecutively sampled    %
%                       around the center of K-space                      %
%           - motion_level: amplitude of rigid fetal movements to         %
%                           simulate                                      %
%           - zip: scanner zero-interpolation filling                     %
%           - ReconMatrix: size of the reconstruction matrix (in voxels)  %
%           - std_noise: standard deviation of the Gaussian noise added   %
%                        to the Fourier domain of the simulated images    %
%           - output_folder: folder where the simulated images are saved  %
%                                                                         %
%  output:  - FSE_Images: simulated T2-weighted MR images of the fetal    %
%                         brain based on the acquisition scheme of FSE    %
%                         sequences                                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2023-02-22                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function FSEimages = FaBiAN_main_CHUV_DA(FetalBrainModelPath, ...
                                                  FetalModel, ...
                                                       SubID, ...
                                                       SesID, ...
                                                       RunID, ...
                                                    Shift_mm, ...
                                                 Orientation, ...
                                                         INU, ...
                                              SamplingFactor, ...
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
                                                   FlipAngle, ...
                                                         ACF, ...
                                                    RefLines, ...
                                                 MotionLevel, ...
                                                         ZIP, ...
                                                 ReconMatrix, ...
                                                     SDnoise, ...
                                               SimResampling, ...
                                                     SimCrop, ...
                                                OutputFolder, ...
                                             WMheterogeneity, ...
                                                       T1_WM, ...
                                                       T2_WM, ...
                                                       T1_GM, ...
                                                       T2_GM, ...
                                                      T1_CSF, ...
                                                      T2_CSF)

% Input check
if nargin < 38
    error('Missing input(s).');
elseif nargin > 38
    error('Too many inputs.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Create output folders / filenames                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create output folders to organize the data into BIDS
% output_path = strcat(output_folder, char(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/');
OutputPath = strcat(OutputFolder, 'data/', num2str(SubID), '/ses-', sprintf('%02s', num2str(SesID)), '/anat/');
OutputPathReo = strcat(OutputFolder, 'data_reo/', num2str(SubID), '/ses-', sprintf('%02s', num2str(SesID)), '/anat/');
if not(isfolder(OutputPath))
    mkdir(OutputPath)
end
if not(isfolder(OutputPathReo))
    mkdir(OutputPathReo)
end
DerivativesPath = strcat(OutputFolder, 'data/derivatives/');
DerivativesPathReo = strcat(OutputFolder, 'data_reo/derivatives/');
if not(isfolder(DerivativesPath))
    mkdir(DerivativesPath)
end
if not(isfolder(DerivativesPathReo))
    mkdir(DerivativesPathReo)
end
% derivatives_labels_path = strcat(derivatives_path, 'labels/', char(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/');
DerivativesLabelsPath = strcat(DerivativesPath, 'labels/', num2str(SubID), '/ses-', sprintf('%02s', num2str(SesID)), '/anat/');
DerivativesLabelsPathReo = strcat(DerivativesPathReo, 'labels/', num2str(SubID), '/ses-', sprintf('%02s', num2str(SesID)), '/anat/');
% derivatives_masks_path = strcat(derivatives_path, 'masks/', char(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/');
DerivativesMasksPath = strcat(DerivativesPath, 'masks/', num2str(SubID), '/ses-', sprintf('%02s', num2str(SesID)), '/anat/');
DerivativesMasksPathReo = strcat(DerivativesPathReo, 'masks/', num2str(SubID), '/ses-', sprintf('%02s', num2str(SesID)), '/anat/');
if not(isfolder(DerivativesLabelsPath))
    mkdir(DerivativesLabelsPath)
end
if not(isfolder(DerivativesLabelsPathReo))
    mkdir(DerivativesLabelsPathReo)
end
if not(isfolder(DerivativesMasksPath))
    mkdir(DerivativesMasksPath)
end
if not(isfolder(DerivativesMasksPathReo))
    mkdir(DerivativesMasksPathReo)
end

% output_im = strcat(output_path, char(sub_id), '_ses-', sprintf('%02s', num2str(ses_id)), '_run-', num2str(run_id), '_T2w.nii');
OutputIm = strcat(OutputPath, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_T2w.nii');
OutputImReo = strcat(OutputPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_T2w_reo.nii');
% OutputImResampled = strcat(OutputPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_T2w_resampled.nii');
OutputImCrop = strcat(OutputPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_T2w_crop.nii');
% Derivatives
OutputLabels = strcat(DerivativesLabelsPath, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_labels.nii');
OutputLabelsReo = strcat(DerivativesLabelsPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_labels_reo.nii');
% OutputLabelsResampled = strcat(DerivativesLabelsPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_labels_resampled.nii');
OutputLabelsCrop = strcat(DerivativesLabelsPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_labels_crop.nii');
OutputMask = strcat(DerivativesMasksPath, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_mask.nii');
OutputMaskReo = strcat(DerivativesMasksPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_mask_reo.nii');
% OutputMaskResampled = strcat(DerivativesMasksPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_mask_resampled.nii');
OutputMaskCrop = strcat(DerivativesMasksPathReo, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_mask_crop.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load fetal brain model and intensity non-uniformity fields             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load segmented high-resolution images of the fetal brain
switch FetalModel
    case 'STA'
        % Gestational age (in weeks)
        GA = SubID;
%         % Session ID
%         ses_id = 1;
    case 'FeTA_CHUV'
        ParticipantsMetadata = strcat(FetalBrainModelPath, 'FeTA-CHUV_participants.xlsx');
        Participants = readtable(ParticipantsMetadata);
        index = string(Participants{:,5})==SubID;
        % Gestational age (in weeks)
        GA = Participants{index,3};
        % Session ID
%         ses_id = participants{index,2};
    case 'FeTA'
        ParticipantsMetadata = strcat(FetalBrainModelPath, '202303-FeTA2021_Release1and2Corrected_v4-Replicate_OHBM_training_set.xlsx');
        Participants = readtable(ParticipantsMetadata, 'Sheet', 'data2sim_1', 'Range', 'A1:F11');
        index = string(Participants{:,1})==SubID;
        % Gestational age (in weeks)
        GA = round(Participants{index,3});
        % Session ID
%         ses_id = unique(participants{index,5});
end

% Load segmented high-resolution anatomical MR images of the fetal brain at
% gestational age GA
[FetalBrain, ModelNiiinfo] = brain_model(FetalBrainModelPath, ...
                                                  FetalModel, ...
                                                       SubID);

% Read the resolution of the 3D anatomical model
SimRes = ModelNiiinfo.PixelDimensions;

% Axis direction codes for the affine orientation matrix of the anatomical
% model of the fetal brain
Affine = ModelNiiinfo.Transform.T';
ModelAxcodes = aff2axcodes(Affine);

% The slice thickness direction will be encoded in the 3. dimension
switch Orientation
    case 1  %sagittal: slice thickness L-R
        SliceDir_idx = find(ismember(ModelAxcodes, ["R","L"]));
        InPlane_idx = find(ismember(ModelAxcodes, ["A","P","S","I"]));
    case 2  %coronal: slice thickness A-P
        SliceDir_idx = find(ismember(ModelAxcodes, ["A","P"]));
        InPlane_idx = find(ismember(ModelAxcodes, ["R","L","S","I"]));
    case 3  %axial: slice thickness S-I
        SliceDir_idx = find(ismember(ModelAxcodes, ["S","I"]));
        InPlane_idx = find(ismember(ModelAxcodes, ["R","L","A","P"]));
end
AxcodesReo = ModelAxcodes([InPlane_idx SliceDir_idx]);

% Reorient the original 3D anatomical model of the fetal brain so that the
% slice thickness direction is encoded in the 3. dimension
[FetalBrainReo, ...
     AffineReo, ...
     SimResReo] = reorient_volume(FetalBrain, ...
                                      Affine, ...
                                      SimRes, ...
                                  AxcodesReo);

% Convert the shift variable in the slice thickness direction into a number
% of voxels
Shift = round(Shift_mm / SimResReo(3));  %in voxels

% Define the slice slab that covers the fetal brain volume after shifting
% the field-of-view in the slice thickness direction
FetalBrainFOV = FOV_shift(FetalBrainReo, Shift, Orientation);

% Update the orientation matrix of the modified fetal brain model
AffineFOV = update_affine(     FetalBrainReo, ...
                                   AffineReo, ...
                               FetalBrainFOV, ...
                          AffineReo(1:3,1:3));

% Load the intensity non-uniformity fields: the third dimension corresponds
% to the slice thickness orientation
B1Map = brainWeb_inu(INU, FetalBrain);

% Apply the same FOV shift to the b1 map as to the fetal brain model
B1MapFOV = FOV_shift(B1Map, Shift, Orientation);

% Reorient the b1 map to match the orientation of the slice slab
[B1Map_reo, ~] = reorient_volume(  B1MapFOV, ...
                                     Affine, ...
                                     SimRes, ...
                                 AxcodesReo);

% Original anatomical models may be up-sampled in the slice thickness
% direction in order to ensure generalizability in the choice of the slice
% thickness, slice gap, etc.
SubunitRes = SimResReo(3) / SamplingFactor;   %mm
% if mod(SliceThickness,SubunitRes)>1e-4
%     error('The subunit resolution (%.1f mm) must be a multiple of the slice thickness, i.e., %.1f mm.\n', SubunitRes, SliceThickness);
% end

FetalBrainUpsampled = sampling_OoP( FetalBrainFOV, ...
                                   SamplingFactor, ...
                                        'nearest');

% Update the orientation matrix of the modified fetal brain model
AffineUpsampled = update_affine(      FetalBrainFOV, ...
                                          AffineFOV, ...
                                FetalBrainUpsampled, ...
                                 AffineFOV(1:3,1:3));

B1MapUpsampled = sampling_OoP(     B1Map_reo, ...
                              SamplingFactor, ...
                                    'linear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Conversion to MR contrast                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Consider only non-zero labels
FetalBrainLabels = unique(FetalBrain);
FetalBrainLabels = FetalBrainLabels(FetalBrainLabels > 0);
% Write labels of the different segmented brain tissues in a list
FetalBrainTissues = permute(FetalBrainLabels, [2,1]);

% Define the interpolation method used to upsample the partial volume maps
% of the WM mask. A bilinear interpolation significantly increases the
% computational burden of the EPG simulations. Therefore, a
% nearest-neighbor interpolation is used to simulate many subjects in a row
InterpolationMethod = 'nearest';
% Generate reference T1 and T2 maps of the fetal brain with the third
% dimension being the slice thickness direction
[RefT1map, RefT2map] = tissue_to_MR(         FetalModel, ...
                                    FetalBrainUpsampled, ...
                                      FetalBrainTissues, ...
                                                     GA, ...
                                                  SubID, ...
                                        WMheterogeneity, ...
                                                 Affine, ...
                                                 SimRes, ...
                                             AxcodesReo, ...
                                                  Shift, ...
                                            Orientation, ...
                                         SamplingFactor, ...
                                    InterpolationMethod, ...
                                                  T1_WM, ...
                                                  T2_WM, ...
                                                  T1_GM, ...
                                                  T2_GM, ...
                                                 T1_CSF, ...
                                                 T2_CSF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extended Phase Graph (EPG) simulations                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run EPG simulations in every voxel of the fetal brain volume
T2decay = compute_t2decay(FetalBrainUpsampled, ...
                               B1MapUpsampled, ...
                                     RefT1map, ...
                                     RefT2map, ...
                                          ETL, ...
                                    FlipAngle, ...
                                          ESP, ...
                               SamplingFactor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parameters of fetal brain FSE acquisition schemes                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate effective FOVPhase based on PhaseOversampling
FOVPhaseOversampling = FOVPhase * (1 + PhaseOversampling);

% Rounding of the number of phase-encoding lines
nPE = round(PhaseResolution * BaseResolution * (1 + PhaseOversampling));
% For Siemens, it probably needs to be rounded to integer values of 2
if mod(nPE,2)~=0
    nPE = nPE + 1;
end

% Calculate a gaussian slice profile as estimated in the SR recon pipeline:
% a Gaussian function with the full width at half maximum equal to the
% slice thickness in the slice-select direction.
SliceProfile = gausswin(round(SliceThickness/SubunitRes), 2*sqrt(2*log(2))/(SliceThickness/SimResReo(3)));

% No slice profile in between slices
Sl_to_Sl = SliceProfile;
Sl_to_Sl(length(SliceProfile)+1:length(SliceProfile)+round(SliceGap/SubunitRes)) = ones(round(SliceGap/SubunitRes),1);

% The number of slices is currently set to ensure that the maximum size in
% mm of the image matrix is an integer number of the slice thickness +
% slice gap
NbSlices = ceil(max([size(T2decay,1), size(T2decay,2), size(T2decay,3)/SamplingFactor]) * SimResReo(3) / SubunitRes / length(Sl_to_Sl));
% Corresponding matrix size
T2decayMaxDim = ceil(NbSlices * length(Sl_to_Sl) / (SimResReo(3) / SubunitRes));

% Zero-padding of the T2decay array so that it is 3D isotropic - This makes
% some calculations below a bit easier
T2decayZP = Resize_Volume(T2decay, [T2decayMaxDim, T2decayMaxDim, T2decayMaxDim*SamplingFactor, size(T2decay,4)]);
% There might be a couple of slices that are just zeros due to this
% zero-padding step.

% In the same way, zero-padding of the segmented fetal brain images after
% upsampling
FetalBrainZP = Resize_Volume(FetalBrainUpsampled, [T2decayMaxDim, T2decayMaxDim, T2decayMaxDim*SamplingFactor]);

% Update the orientation matrix of the modified fetal brain model while
% taking into account the deviation from the center of the anatomical model
% of the fetal brain which may result from this zero-padding step
CenterOffset = [0 0 0];
if mod(size(FetalBrainZP,1)-size(FetalBrainUpsampled,1),2)~=0
    if size(FetalBrainZP,1)>size(FetalBrainUpsampled,1)
        CenterOffset(1) = 1;    %The center of the volume is updated by '+1' voxel in this direction because we zero-padd the volume, meaning the final size is larger than the original one
    else
        CenterOffset(1) = -1;
    end
end
if mod(size(FetalBrainZP,2)-size(FetalBrainUpsampled,2),2)~=0
    if size(FetalBrainZP,2)>size(FetalBrainUpsampled,2)
        CenterOffset(2) = 1;
    else
        CenterOffset(2) = -1;
    end
end
if mod(size(FetalBrainZP,3)-size(FetalBrainUpsampled,3),2)~=0
    if size(FetalBrainZP,3)>size(FetalBrainUpsampled,3)
        CenterOffset(3) = 1;
    else
        CenterOffset(3) = -1;
    end
end
AffineZP = update_affine(         FetalBrainUpsampled, ...
                                      AffineUpsampled, ...
                                         FetalBrainZP, ...
                             AffineUpsampled(1:3,1:3), ...
                         'CenterOffset', CenterOffset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  K-space sampling of the simulated FSE images                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data sampling follows an interleaved acquisition scheme
InterleavedSlices_index = interleaved_scheme(NbSlices);

[  MotionCorruptedSlices, ...
 TranslationDisplacement, ...
           RotationAngle, ...
            RotationAxis] = motion_transform(MotionLevel, ...
                                                NbSlices);

[KSpace, SlVolume] = Kspace_sampling(              T2decayZP, ...
                                               FetalBrainFOV, ...
                                                  FetalModel, ...
                                              B1MapUpsampled, ...
                                           FetalBrainTissues, ...
                                                   SimResReo, ...
                                              SamplingFactor, ...
                                       MotionCorruptedSlices, ...
                                     TranslationDisplacement, ...
                                               RotationAngle, ...
                                                RotationAxis, ...
                                                   'nearest', ...
                                                         ETL, ...
                                                   FlipAngle, ...
                                                         ESP, ...
                                                       TEeff, ...
                                                          TR, ...
                                                    NbSlices, ...
                                     InterleavedSlices_index, ...
                                                    Sl_to_Sl, ...
                                                SliceProfile, ...
                                              SliceThickness, ...
                                                     FOVRead, ...
                                        FOVPhaseOversampling, ...
                                              BaseResolution, ...
                                                         nPE, ...
                                                         ACF, ...
                                                    RefLines, ...
                                                          GA, ...
                                                       SubID, ...
                                             WMheterogeneity, ...
                                                      Affine, ...
                                                      SimRes, ...
                                                  AxcodesReo, ...
                                                       Shift, ...
                                                 Orientation, ...
                                                       T1_WM, ...
                                                       T2_WM, ...
                                                       T1_GM, ...
                                                       T2_GM, ...
                                                      T1_CSF, ...
                                                      T2_CSF);

SqueezeSlVolume = squeeze(SlVolume(:,:,:,30));
SlVolume_resized = Resize_Volume(SqueezeSlVolume, [round(FOVRead/SimResReo(1)), round(FOVPhaseOversampling/SimResReo(2)), NbSlices]);

% Update the orientation matrix of the modified fetal brain model while
% taking into account the deviation from the center of the anatomical model
% of the fetal brain which may result from this zero-padding step
CenterOffset = [0 0 0];
if mod(size(SlVolume_resized,1)-size(SqueezeSlVolume,1),2)~=0
    if size(SlVolume_resized,1)>size(SqueezeSlVolume,1)
        CenterOffset(1) = 1;
    else
        CenterOffset(1) = -1;
    end
end
if mod(size(SlVolume_resized,2)-size(SqueezeSlVolume,2),2)~=0
    if size(SlVolume_resized,2)>size(SqueezeSlVolume,2)
        CenterOffset(2) = 1;
    else
        CenterOffset(2) = -1;
    end
end
if mod(size(SlVolume_resized,3)-size(SqueezeSlVolume,3),2)~=0
    if size(SlVolume_resized,3)>size(SqueezeSlVolume,3)
        CenterOffset(3) = 1;
    else
        CenterOffset(3) = -1;
    end
end
% Voxel size of the simulated images
SimVoxSize = [FOVRead/ReconMatrix (FOVPhaseOversampling/(1+PhaseOversampling))/ReconMatrix SliceThickness+SliceGap];
% Initialize the rotation matrix
RotationsSlVolume = (AffineZP(1:3, 1:3) ./ SimResReo) .* SimVoxSize;
AffineSlVolume = update_affine(                FetalBrainZP, ...
                                                   AffineZP, ...
                                           SlVolume_resized, ...
                                          RotationsSlVolume, ...
                               'CenterOffset', CenterOffset);

% Automatically propagate labels of the simulated images
SimLabels = auto_segmentation(          FetalBrainZP, ...
                                       FetalBrainReo, ...
                                           SimResReo, ...
                                      SamplingFactor, ...
                                      SliceThickness, ...
                                          SubunitRes, ...
                                             FOVRead, ...
                                FOVPhaseOversampling, ...
                                            NbSlices, ...
                             InterleavedSlices_index, ...
                                            Sl_to_Sl, ...
                               MotionCorruptedSlices, ...
                             TranslationDisplacement, ...
                                       RotationAngle, ...
                                        RotationAxis);

% Reduce the labels from Gholipour atlas to 6 main labels
switch FetalModel
    case 'STA'
    % Group structures and assign them the same labels as in the FeTA dataset +
    % merge extra-axial CSF space (class 1) and intra-cranial CSF, ventricles
    % system (class 4) into one single class (class 1)
    %     0 Background: [0]
    %     1 Extra-axial CSF space: [124]
    %     2 Gray matter: [37,38,41,42,112,113]
    % NB-FeTA annotation guidelines: "The label includes the amygdala and
    % hippocampus, but does not include any infratentorial structures."
    %     3 White matter: [91,110,111,114,115,116,117,118,119,120,121,125]
    %     4 Ventricular system: [92,93]
    %     5 Cerebellum: [100,101]
    %     6 Deep gray matter: [71,72,73,74,77,78,108,109,122,123]
    % NB-FeTA annotation guidelines: "The internal capsule is within this
    % label, as it is often not possible to separate due to low contrast in the
    % 3D reconstructed T2 images.
    %     7 Brainstem: [94]
    SimLabels(SimLabels(:)==124) = 1;
    SimLabels(SimLabels(:)==37|SimLabels(:)==38|SimLabels(:)==41|SimLabels(:)==42|SimLabels(:)==112|SimLabels(:)==113) = 2;
    SimLabels(SimLabels(:)==91|SimLabels(:)==110|SimLabels(:)==111|SimLabels(:)==114|SimLabels(:)==115|SimLabels(:)==116|SimLabels(:)==117|SimLabels(:)==118|SimLabels(:)==119|SimLabels(:)==120|SimLabels(:)==121|SimLabels(:)==125) = 3;
    SimLabels(SimLabels(:)==92|SimLabels(:)==93) = 1;
    SimLabels(SimLabels(:)==100|SimLabels(:)==101) = 5;
    SimLabels(SimLabels(:)==71|SimLabels(:)==72|SimLabels(:)==73|SimLabels(:)==74|SimLabels(:)==77|SimLabels(:)==78|SimLabels(:)==108|SimLabels(:)==109|SimLabels(:)==122|SimLabels(:)==123) = 6;
    SimLabels(SimLabels(:)==94) = 7;
end

% Simulate scanner zero-interpolation filling (ZIP)
if ZIP==1   %ZIP
    KSpace = zip_kspace(                                                                                              KSpace, ...
                        [size(KSpace,1)*ReconMatrix/BaseResolution size(KSpace,2)*ReconMatrix/BaseResolution size(KSpace,3)]);
end

% Add complex Gaussian noise to K-space
KSpaceNoise = add_noise(KSpace, SDnoise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Final simulated FSE images                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn back data in K-space to the image space
ComplexFSEimages = imresize(ifft2c(KSpaceNoise), [size(KSpace,1), round(size(KSpace,2)/PhaseResolution)]);

% Histogram normalization between 0 and 255 of simulated SS-FSE images
NormFSEimages = ComplexFSEimages * 255 / max(abs(ComplexFSEimages(:)));

% Only keep the magnitude image
FSEimages = abs(NormFSEimages);

% Associated label map
LabelMap = imresize(SimLabels, [size(KSpace,1), round(size(KSpace,2)/PhaseResolution)], 'nearest');

% Resize the simulated images to the desired dimensions (especially, do not
% reconstruct the phase oversampling)
AcqDim = [ReconMatrix ReconMatrix NbSlices];
if any(size(FSEimages))~=AcqDim
    ReconFSEimages = Resize_Volume(FSEimages, AcqDim);
    ReconLabelMap = Resize_Volume(LabelMap, AcqDim);
else
    ReconFSEimages = FSEimages;
    ReconLabelMap = LabelMap;
end

% Initialize the orientation matrix of the simulated images
% NewRotations = (Affine_SlVolume(1:3, 1:3) ./ SimResReo) .* SimVoxSize;
% Update the affine matrix of the modified fetal brain model while taking
% into account any deviation from the center of the anatomical model of the
% fetal brain which may result from zero-padding when sampling K-space
CenterOffset = [0 0 0];
if mod(AcqDim(1)-size(SlVolume_resized,1),2)~=0
    if AcqDim(1)>size(SlVolume_resized,1)
        CenterOffset(1) = 1;
    else
        CenterOffset(1) = -1;
    end
end
if mod(AcqDim(2)-size(SlVolume_resized,2),2)~=0
    if AcqDim(2)>size(SlVolume_resized,2)
        CenterOffset(2) = 1;
    else
        CenterOffset(2) = -1;
    end
end
if mod(AcqDim(3)-size(SlVolume_resized,3),2)~=0
    if AcqDim(3)>size(SlVolume_resized,3)
        CenterOffset(3) = 1;
    else
        CenterOffset(3) = -1;
    end
end
ReconAffine = update_affine(            SlVolume_resized, ...
                                          AffineSlVolume, ...
                                          ReconFSEimages, ...
                                 AffineSlVolume(1:3,1:3), ...
                            'CenterOffset', CenterOffset);

% Write nifti header of the simulated images with the slice thickness
% encoded in the third dimension
ReconNiiinfo = ModelNiiinfo;
ReconNiiinfo.Transform.T = ReconAffine';
ReconNiiinfo.TransformName = "Sform";
ReconNiiinfo.ImageSize = size(ReconFSEimages);
ReconNiiinfo.PixelDimensions = [sqrt(sum(ReconAffine(1:3,1).^2)) sqrt(sum(ReconAffine(1:3,2).^2)) sqrt(sum(ReconAffine(1:3,3).^2))];
if isa(ReconFSEimages, ModelNiiinfo.Datatype)==0
    ReconNiiinfo.Datatype = class(ReconFSEimages);
end

% Save the simulated images after reorientation (i.e., the slice thickness
% is encoded in the third dimension as for clinical acquisitions) in nifti
% format
niftiwrite(ReconFSEimages, OutputImReo, ReconNiiinfo, 'Compressed', true);

% Generate a binary mask from the reconstructed label map
ReconBinaryMask = ReconLabelMap;
ReconBinaryMask(ReconBinaryMask(:)~=0) = 1;

% Save corresponding binary masks and label maps after reorientation (i.e.,
% the slice thickness is encoded in the third dimension as for clinical
% acquisitions) in nifti format
MaskNiiinfo = ReconNiiinfo;
if isa(ReconLabelMap, ReconNiiinfo.Datatype)==0
    MaskNiiinfo.Datatype = class(ReconLabelMap);
end
niftiwrite(ReconLabelMap, OutputLabelsReo, MaskNiiinfo, 'Compressed', true);
niftiwrite(ReconBinaryMask, OutputMaskReo, MaskNiiinfo, 'Compressed', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Depending on the application, resampling of the simulated images       %
%  might be needed                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch SimResampling
    case 'true'
        % Downsample the simulated SS-FSE images, the corresponding label map and
        % mask in the in-plane direction to match the in-plane resolution of the
        % clinical SR reconstructions
        ResamplingRead = SimResReo(1) / (FOVRead/ReconMatrix);
        ResamplingPhase = SimResReo(2) / (FOVPhase/ReconMatrix);
        % Resample simulated images and derivatives
        ResReconImages = resampling_inplane( ReconFSEimages, ...
                                             ResamplingRead, ...
                                            ResamplingPhase, ...
                                                   'linear');
        ResReconLabels = resampling_inplane(  ReconLabelMap, ...
                                             ResamplingRead, ...
                                            ResamplingPhase, ...
                                                  'nearest');
        % Generate a binary mask from the reconstructed label map
        ResReconMask = ResReconLabels;
        ResReconMask(ResReconMask(:)~=0) = 1;
        % Voxel size of the reconstructed simulated images
        SimVoxSize = [(FOVRead/ReconMatrix)*ResamplingRead (FOVPhaseOversampling/(1+PhaseOversampling))/ReconMatrix*ResamplingPhase SliceThickness+SliceGap];
        % Initialize the orientation matrix of the simulated images
        NewRotations = ReconAffine(1:3,1:3) .* [ResamplingRead ResamplingPhase 1];
        % Update the orientation matrix of the modified fetal brain model
        CenterOffset = [0 0 0];
        if mod(size(ResReconImages,1)-size(ReconFSEimages,1),2)~=0
            if size(ResReconImages,1)>size(ReconFSEimages,1)
                CenterOffset(1) = 1;
            else
                CenterOffset(1) = -1;
            end
        end
        if mod(size(ResReconImages,2)-size(ReconFSEimages,2),2)~=0
            if size(ResReconImages,2)>size(ReconFSEimages,2)
                CenterOffset(2) = 1;
            else
                CenterOffset(2) = -1;
            end
        end
        if mod(size(ResReconImages,3)-size(ReconFSEimages,3),2)~=0
            if size(ResReconImages,3)>size(ReconFSEimages,3)
                CenterOffset(3) = 1;
            else
                CenterOffset(3) = -1;
            end
        end
        ResReconAffine = update_affine(              ReconFSEimages, ...
                                                        ReconAffine, ...
                                                     ResReconImages, ...
                                                       NewRotations, ...
                                       'CenterOffset', CenterOffset);
        % Save the resampled, reconstructed simulated images after
        % reorientation (i.e., the slice thickness is encoded in the third
        % dimension as for clinical acquisitions) in nifti format
        ResReconNiiinfo = ModelNiiinfo;
        ResReconNiiinfo.Transform.T = ResReconAffine';
        ResReconNiiinfo.TransformName = "Sform";
        ResReconNiiinfo.ImageSize = size(ResReconImages);
        ResReconNiiinfo.PixelDimensions = [sqrt(sum(ResReconAffine(1:3,1).^2)) sqrt(sum(ResReconAffine(1:3,2).^2)) sqrt(sum(ResReconAffine(1:3,3).^2))];
        if isa(ResReconImages, ModelNiiinfo.Datatype)==0
            ResReconNiiinfo.Datatype = class(ResReconImages);
        end
        niftiwrite(ResReconImages, OutputImResampled, ResReconNiiinfo, 'Compressed', true);
        ResMaskNiiinfo = ResReconNiiinfo;
        if isa(ResReconLabels, ResReconNiiinfo.Datatype)==0
            ResMaskNiiinfo.Datatype = class(ResReconLabels);
        end
        niftiwrite(ResReconLabels, OutputLabelsResampled, ResMaskNiiinfo, 'Compressed', true);
        niftiwrite(ResReconMask, OutputMaskResampled, ResMaskNiiinfo, 'Compressed', true);
    case 'false'
%         % Voxel size of the simulated images
%         SimVoxSize = [FOVRead/ReconMatrix (FOVPhase_oversampling/(1+PhaseOversampling))/ReconMatrix SliceThickness+SliceGap];
        ResReconAffine = ReconAffine;
        ResReconImages = ReconFSEimages;
        ResReconLabels = ReconLabelMap;
        ResReconMask = ReconBinaryMask;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Save the final simulated FSE images in the same orientation and space  %
%  coordinates as the original anatomical model                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Come back to the original fetal brain model space coordinates
[SimFSEImages, SimAffineModelAxcodes, ~] = reorient_volume(ResReconImages, ...
                                                           ResReconAffine, ...
                                                                SimResReo, ...
                                                             ModelAxcodes);
[SimLabelMap, ~] = reorient_volume(ResReconLabels, ...
                                      ReconAffine, ...
                                        SimResReo, ...
                                     ModelAxcodes);
[SimBinaryMask, ~] = reorient_volume( ResReconMask, ...
                                       ReconAffine, ...
                                         SimResReo, ...
                                      ModelAxcodes);
% sim_axcodes = aff2axcodes(affine_reo);

% Write nifti header of the simulated images after reorientation in the
% original coordinate system of the anatomical model
SimModelNiiinfo = ModelNiiinfo;
SimModelNiiinfo.Transform.T = SimAffineModelAxcodes';
% SimModelNiiinfo.TransformName = "Qform";
SimModelNiiinfo.TransformName = "Sform";
SimModelNiiinfo.ImageSize = size(SimFSEImages);
SimModelNiiinfo.PixelDimensions = [sqrt(sum(SimAffineModelAxcodes(1:3, 1).^2)) sqrt(sum(SimAffineModelAxcodes(1:3, 2).^2)) sqrt(sum(SimAffineModelAxcodes(1:3, 3).^2))];
if isa(SimFSEImages, ModelNiiinfo.Datatype)==0
    SimModelNiiinfo.Datatype = class(SimFSEImages);
end

% Save simulated images in nifti format, in the same orientation and space
% coordinates as the original anatomical model
niftiwrite(SimFSEImages, OutputIm, SimModelNiiinfo, 'Compressed', true);
% Save corresponding binary masks and label maps in nifti format, in the
% same orientation and space coordinates as the original anatomical model
SimMaskNiiinfo = SimModelNiiinfo;
if isa(SimLabelMap, SimModelNiiinfo.Datatype)==0
    SimMaskNiiinfo.Datatype = class(SimLabelMap);
end
niftiwrite(SimLabelMap, OutputLabels, SimMaskNiiinfo, 'Compressed', true);
niftiwrite(SimBinaryMask, OutputMask, SimMaskNiiinfo, 'Compressed', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Optional: Crop the simulated LR series to the same size as the         %
%  original anatomical model (after reorientation, i.e., the slice        %
%  thickness is encoded in the 3. dimension) and save them as nifti       %
%  after updating the nifti header information                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SimCrop=="true"
    % Crop the simulated images and corresponding derivatives
    CropFSEimages = Resize_Volume(ResReconImages, size(FetalBrainReo));
    CropLabels = Resize_Volume(ResReconLabels, size(FetalBrainReo));
    CropMask = Resize_Volume(ResReconMask, size(FetalBrainReo));
    % Update the orientation matrix of the modified fetal brain model to
    % account for this additional zero-padding step
    CenterOffset = [0 0 0];
    if mod(size(CropFSEimages,1)-size(ResReconImages,1),2)~=0
        if size(CropFSEimages,1)>size(ResReconImages,1)
            CenterOffset(1) = 1;
        else
            CenterOffset(1) = -1;   %Here, the center should be deviated by '-1' because we crop the image (i.e., the size of the new volume is smaller than the original one)
        end
    end
    if mod(size(CropFSEimages,2)-size(ResReconImages,2),2)~=0
        if size(CropFSEimages,2)>size(ResReconImages,2)
            CenterOffset(2) = 1;
        else
            CenterOffset(2) = -1;
        end
    end
    if mod(size(CropFSEimages,3)-size(ResReconImages,3),2)~=0
        if size(CropFSEimages,3)>size(ResReconImages,3)
            CenterOffset(3) = 1;
        else
            CenterOffset(3) = -1;
        end
    end
    CropAffine = update_affine(              ResReconImages, ...
                                             ResReconAffine, ...
                                              CropFSEimages, ...
                                    ResReconAffine(1:3,1:3), ...
                               'CenterOffset', CenterOffset);
    % Write nifti header of the cropped simulated images with the slice
    % thickness encoded in the third dimension
    ModelCropNiiinfo = ModelNiiinfo;
    ModelCropNiiinfo.Transform.T = CropAffine';
    ModelCropNiiinfo.TransformName = "Sform";
    ModelCropNiiinfo.ImageSize = size(CropFSEimages);
    ModelCropNiiinfo.PixelDimensions = [sqrt(sum(CropAffine(1:3, 1).^2)) sqrt(sum(CropAffine(1:3, 2).^2)) sqrt(sum(CropAffine(1:3, 3).^2))];
    if isa(CropFSEimages, ModelNiiinfo.Datatype)==0
        ModelCropNiiinfo.Datatype = class(CropFSEimages);
    end
    % Save the simulated cropped images after reorientation (i.e., the
    % slice thickness is encoded in the third dimension as for clinical
    % acquisitions) in nifti format
    niftiwrite(CropFSEimages, OutputImCrop, ModelCropNiiinfo, 'Compressed', true);
    % Save the corresponding binary mask and label map in nifti format
    CropMaskNiiinfo = ModelCropNiiinfo;
    if isa(CropLabels, ModelCropNiiinfo.Datatype)==0
        CropMaskNiiinfo.Datatype = class(CropLabels);
    end
    niftiwrite(CropLabels, OutputLabelsCrop, CropMaskNiiinfo, 'Compressed', true);
    niftiwrite(CropMask, OutputMaskCrop, CropMaskNiiinfo, 'Compressed', true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Save raw data                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save(strcat(output_folder, 'KSpace'), 'KSpace', '-v7.3');
% save(strcat(output_folder, 'KSpace_noise'), 'KSpace_noise', '-v7.3');
% save(strcat(output_folder, 'Fetal_Brain_FSE_Images'), 'FSE_Images', '-v7.3');
% save(strcat(output_folder, 'Label_Map_recon'), 'Label_Map_recon', '-v7.3');
% save(strcat(output_folder, 'Motion_Transforms'), 'motion_corrupted_slices', 'translation_displacement', 'rotation_angle', 'rotation_axis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Save the parameters of the MR acquisition simulated                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save simulation parameters
save('sim_parameters.mat',                'SubID', ...
                                          'SesID', ...
                                          'RunID', ...
                                     'FetalModel', ...
                                    'Orientation', ...
                                       'Shift_mm', ...
                                            'INU', ...
                                 'SamplingFactor', ...
                                             'B0', ...
                                            'ESP', ...
                                            'ETL', ...
                              'PhaseOversampling', ...
                                 'SliceThickness', ...
                                       'SliceGap', ...
                                        'FOVRead', ...
                                       'FOVPhase', ...
                           'FOVPhaseOversampling', ...
                                 'BaseResolution', ...
                                'PhaseResolution', ...
                                             'TR', ...
                                          'TEeff', ...
                                      'FlipAngle', ...
                                            'ACF', ...
                                       'RefLines', ...
                                    'MotionLevel', ...
                                            'ZIP', ...
                                    'ReconMatrix', ...
                                  'SimResampling', ...
                                     'SimVoxSize', ...
                                        'SDnoise', ...
                                'WMheterogeneity')

if SimResampling=="true"
    save('resampling.mat', 'ResamplingRead', 'ResamplingPhase')
end

% Write simulations parameters in a .csv file
data = load('sim_parameters.mat');
param = fieldnames(data);
% writecell(horzcat(param, struct2cell(data)), strcat(output_folder, char(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/', char(sub_id), '_ses-', sprintf('%02s', num2str(ses_id)), '_run-', num2str(run_id), '_sim_parameters.csv'))
writecell(horzcat(param, struct2cell(data)), strcat(OutputPath, num2str(SubID), '_ses-', sprintf('%02s', num2str(SesID)), '_run-', num2str(RunID), '_sim_parameters.csv'))
% writecell(horzcat(param, struct2cell(data)), strcat(output_folder, 'data_reo/', num2str(sub_id), '/ses-', sprintf('%02s', num2str(ses_id)), '/anat/', num2str(sub_id), '_ses-', sprintf('%02s', num2str(ses_id)), '_run-', num2str(run_id), '_sim_parameters.csv'))

end