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
%                                                 sampling_factor, ...    %
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
%                                                     reconMatrix, ...    %
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
%           - sampling_factor: factor of up- or down-sampling             %
%           - B0: main magnetic field strength                            %
%           - ESP: echo spacing (in ms)                                   %
%           - ETL: echo train length, i.e. number of 180°-RF pulses       %
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
%           - ACF: acceleration factor for K-space sampling,              %
%                 in Half Acquisition Fourier is 2                        %
%           - RefLines: number of lines that are consecutively sampled    %
%                       around the center of K-space                      %
%           - motion_level: amplitude of rigid fetal movements to         %
%                           simulate                                      %
%           - zip: scanner zero-interpolation filling                     %
%           - reconMatrix: size of the reconstruction matrix (in voxels)  %
%           - std_noise: standard deviation of the Gaussian noise added   %
%                        to the Fourier domain of the simulated images    %
%           - output_folder: folder where the simulated images are saved  %
%                                                                         %
%  output:  - FSE_Images: simulated T2-weighted MR images of the fetal    %
%                         brain based on the acquisition scheme of FSE    %
%                         sequences                                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-22                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function FSE_Images = FaBiAN_main(Fetal_Brain_model_path, ...
                                                      GA, ...
                                                  SimRes, ...
                                                  acq_voxel_size,...
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
                                           output_folder)

% Input check
if nargin < 27
    error('Missing input(s).');
elseif nargin > 27
    error('Too many inputs.');
end

addpath('Utilities')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load fetal brain model and intensity non-uniformity fields             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd .. 
% Load segmented high-resolution images of the fetal brain at gestational
% age GA
Fetal_Brain = brain_model(Fetal_Brain_model_path, GA);

% Convert the shift variable into a number of voxels
shift = shift_mm / SimRes;  %in voxels

% Define the slice slab that covers the fetal brain volume based on the
% shift variable in the slice thickness direction and the acquisition plane
Fetal_Brain = FOV_shift(Fetal_Brain, shift, orientation);

% Load the intensity non-uniformity fields
b1map = brainWeb_inu(inu, Fetal_Brain);

% Reorient volumes
% For the sake of clarity, the orientation corresponding to the slice
% thickness will always be given by the 3. index
Fetal_Brain = volume_reorient(Fetal_Brain, orientation);
b1map = volume_reorient(b1map, orientation);

% Original atlas images are isotropic (SimRes in mm). They are up-sampled
% to 0.1 mm in what will be defined as the slice thickness direction in
% order to be able to accurately account for the slice profile, the slice
% thickness, the slice gap, etc, ..., whatever the situation, keeping the
% framework as general as possible.
SubunitRes = SimRes / sampling_factor;   %mm
resample_size = SimRes./acq_voxel_size;
full_resized= [resample_size(1),resample_size(2),sampling_factor];
Fetal_Brain_upsampled = imresize3(Fetal_Brain, ...
                                     round(size(Fetal_Brain).*full_resized), ...
                                           'bilinear');

b1map_upsampled = imresize3(          b1map, ...
                               round(size(b1map).*full_resized), ...
                                      'bilinear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Conversion to MR contrast                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Consider only non-zero labels
Fetal_Brain_Labels = unique(Fetal_Brain);
Fetal_Brain_Labels = Fetal_Brain_Labels(Fetal_Brain_Labels > 0);
% Write labels of the different segmented brain tissues in a list
Fetal_Brain_Tissues = permute(Fetal_Brain_Labels, [2,1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parameters of fetal brain FSE acquisition schemes                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate effective FOVPhase based on PhaseOversampling
FOVPhase = FOVPhase * (1 + PhaseOversampling);

% Rounding of the number of phase-encoding lines
nPE = round(PhaseResolution * BaseResolution * (1 + PhaseOversampling));
% For Siemens, it probably needs to be rounded to integer values of 2
if mod(nPE,2)~=0
    nPE = nPE + 1;
end

% Calculate a gaussian slice profile as estimated in the SR recon pipeline:
% a Gaussian function with the full width at half maximum equal to the
% slice thickness in the slice-select direction.
SliceProfile = gausswin(SliceThickness/SubunitRes, 2*sqrt(2*log(2))/(SliceThickness/SimRes));

% No slice profile in between slices
Sl_to_Sl = SliceProfile;
Sl_to_Sl(length(SliceProfile)+1:length(SliceProfile)+SliceGap/SubunitRes) = ones(round(SliceGap/SubunitRes),1);

% The number of slices is currently set to ensure that the maximum size in
% mm of the image matrix is an integer number of the slice thickness +
% slice gap
NbSlices = ceil(max([size(Fetal_Brain_upsampled,1), size(Fetal_Brain_upsampled,2), size(Fetal_Brain_upsampled,3)/sampling_factor]) * SimRes / SubunitRes / length(Sl_to_Sl));
% Corresponding matrix size
T2decay_MaxDim = ceil(NbSlices * length(Sl_to_Sl) / (SimRes / SubunitRes));
Fetal_Brain_upsampled=Resize_Volume(Fetal_Brain_upsampled,[T2decay_MaxDim,T2decay_MaxDim,T2decay_MaxDim*sampling_factor]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Motion generation                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data sampling follows an interleaved acquisition scheme
interleavedSlices_index = interleaved_scheme(NbSlices);

[ motion_corrupted_slices, ...
 translation_displacement, ...
           rotation_angle, ...
            rotation_axis] = motion_transform(motion_level, ...
                                                  NbSlices);
Fetal_Brain_rotated=Fetal_Brain;
disp("Simulating brain movement")
motion_index = 0;
for iSlice=1:length(interleavedSlices_index) 
    % Inter-slice motion -
    % - Random translation: uniform distribution between
    % [-translation_amplitude,+translation_amplitude]mm
    if any(interleavedSlices_index(iSlice)==motion_corrupted_slices)
        motion_index = motion_index + 1;
        Fetal_Brain_rotated = moving_3Dbrain(                             Fetal_Brain_rotated, ...
                                                                               SimRes, ...
                                             translation_displacement(motion_index,:), ...
                                                                            'nearest', ...
                                                         rotation_angle(motion_index), ...
                                                        rotation_axis(motion_index,:), ...
                                                                            'nearest', ...
                                                                               'crop');
        Fetal_Brain_rotated_upsampled = sampling_OoP(            Fetal_Brain_rotated, ...
                                                                     sampling_factor, ...
                                                     "linear");%interpolation_method_upsampling;
        Fetal_Brain_rotated_upsampled=Resize_Volume(Fetal_Brain_rotated_upsampled,[T2decay_MaxDim,T2decay_MaxDim,T2decay_MaxDim*sampling_factor]);

        
        sliceIndex=interleavedSlices_index(iSlice)*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
        Fetal_Brain_upsampled(:,:,sliceIndex:sliceIndex+SliceThickness/SubunitRes-1) = Fetal_Brain_rotated_upsampled(:,:,sliceIndex:sliceIndex+SliceThickness/SubunitRes-1);
        %Perform the rotation for all following slices
        if iSlice < length(interleavedSlices_index)
            for nextSlice=iSlice+1:length(interleavedSlices_index)
                nextSlice_index = interleavedSlices_index(nextSlice)*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
                Fetal_Brain_upsampled(:,:,nextSlice_index:nextSlice_index+SliceThickness/SubunitRes-1) = Fetal_Brain_rotated_upsampled(:,:,nextSlice_index:nextSlice_index+SliceThickness/SubunitRes-1);
            end
        end
    end
end
% Updating T1 and T2 maps
[ref_T1map, ref_T2map] = tissue_to_MR(Fetal_Brain_upsampled,Fetal_Brain_Tissues, B0);
b1map_upsampled=Resize_Volume(b1map_upsampled,[T2decay_MaxDim,T2decay_MaxDim,T2decay_MaxDim*sampling_factor]);


% We start the slice-by-slice computation 
KSpace = zeros(BaseResolution, nPE,NbSlices, 'single');
h = waitbar(0,'Computing slice-by-slice KSpace');
counter=0;
for iSlice=interleavedSlices_index
    counter=counter+1;
    index = iSlice*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Extended Phase Graph (EPG) simulations                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nextIndex=index+SliceThickness/SubunitRes-1;

    % Run EPG simulations in every voxel of the fetal brain volume
    T2decay_slice = compute_t2decay(Fetal_Brain_upsampled(:,:,index:nextIndex), ...
                                    b1map_upsampled(:,:,index:nextIndex), ...
                                          ref_T1map(:,:,index:nextIndex), ...
                                          ref_T2map(:,:,index:nextIndex), ...
                                                ETL, ...
                                                ESP, ...
                                    sampling_factor);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  K-space sampling of the simulated FSE images                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data sampling follows an interleaved acquisition scheme


KSpace(:,:,iSlice) = Kspace_sampling(T2decay_slice, ...
                                           SimRes, ...
                                           sampling_factor,...
                                            TEeff, ...
                                               TR, ...
                          iSlice, ...
                                     SliceProfile, ...
                                   SliceThickness, ...
                                          FOVRead, ...
                                         FOVPhase, ...
                                   BaseResolution, ...
                                              nPE, ...
                                              ACF, ...
                                         RefLines);
waitbar(counter/length(interleavedSlices_index))
end
close(h)
% Simulate scanner zero-interpolation filling (ZIP)
if zip==1   %ZIP
    KSpace = zip_kspace(                                  KSpace, ...
                        [reconMatrix reconMatrix size(KSpace,3)]);
end

% Add complex Gaussian noise to K-space
KSpace_noise = add_noise(KSpace, std_noise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Final simulated FSE images                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn back data in K-space to the image space
FSE_Images = imresize(real(ifft2c(KSpace_noise)), [size(KSpace,1), round(size(KSpace,2)/PhaseResolution)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Save data                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save(strcat(output_folder, 'KSpace'), 'KSpace', '-v7.3');
save(strcat(output_folder, 'KSpace_noise'), 'KSpace_noise', '-v7.3');
save(strcat(output_folder, 'Fetal_Brain_FSE_Images'), 'FSE_Images', '-v7.3');
save(strcat(output_folder, 'Motion_Transforms'), 'motion_corrupted_slices', 'translation_displacement', 'rotation_angle', 'rotation_axis');

end