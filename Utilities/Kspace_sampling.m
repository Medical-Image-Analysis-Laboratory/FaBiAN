%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that samples the Fourier domain (K-space) of the simulated    %
%  images according to the acquisition scheme implemented in Fast Spin    %
%  Echo (FSE) sequences.                                                  %
%                                                                         %
%    [   KSpace, ...                                                      %
%     Sl_Volume] = Kspace_sampling(                     T2decay_zp, ...   %
%                                                      Fetal_Brain, ...   %
%                                                      fetal_model, ...   %
%                                                  b1map_upsampled, ...   %
%                                              Fetal_Brain_Tissues, ...   %
%                                                        SimResReo, ...   %
%                                                   SamplingFactor, ...   %
%                                          motion_corrupted_slices, ...   %
%                                         translation_displacement, ...   %
%                                                   rotation_angle, ...   %
%                                                    rotation_axis, ...   %
%                                  interpolation_method_upsampling, ...   %
%                                                              ETL, ...   %
%                                                        flipAngle, ...   %
%                                                              ESP, ...   %
%                                                            TEeff, ...   %
%                                                               TR, ...   %
%                                                         NbSlices, ...   %
%                                          interleavedSlices_index, ...   %
%                                                         Sl_to_Sl, ...   %
%                                                     SliceProfile, ...   %
%                                                   SliceThickness, ...   %
%                                                          FOVRead, ...   %
%                                                         FOVPhase, ...   %
%                                                   BaseResolution, ...   %
%                                                              nPE, ...   %
%                                                              ACF, ...   %
%                                                         RefLines, ...   %
%                                                               GA, ...   %
%                                                           sub_id, ...   %
%                                                 WM_heterogeneity, ...   %
%                                                           affine, ...   %
%                                                           SimRes, ...   %
%                                                          axcodes)       %
%                                                                         %
%  inputs:  - T2decay_zp: 4D matrix that combines the anatomical          %
%                         information from the segmented high-resolution  %
%                         fetal brain volume and the T2 decay computed    %
%                         in every voxel of this volume after zero-       %
%                         -padding                                        %
%           - Fetal_Brain: segmented high-resolution 3D volume of the     %
%                          fetal brain                                    %
%           - fetal_model: anatomical model of the fetal brain            %
%           - b1map_upsampled: B1+ bias field map resampled by the        %
%                              sampling_factor in the slice thickness     %
%                              direction                                  %
%           - Fetal_Brain_Tissues: list of tissues in the fetal brain     %
%                                  that were segmented and labeled        %
%           - SimResReo: resolution of the original 3D high-resolution    %
%                        anatomical model of the fetal brain after        %
%                        reorientation (i.e., the slice thickness is      %
%                        encoded in the third dimension)                  %
%           - SamplingFactor: factor of up- or down-sampling              %
%           - motion_corrupted_slices: list of indexes of the slices      %
%                                      corrupted by motion                %
%           - translation_displacement: table of translation              %
%                                       displacements (in mm) along the   %
%                                       three main axes for each motion-  %
%                                       -corrupted slice                  %
%           - rotation_angle: list of rotation angles (in �) for 3D       %
%                             rotation in every motion-corrupted slice    %
%           - rotation_axis: table of rotation axes for 3D rotation in    %
%                            every motion-corrupted slice                 %
%           - interpolation_method_upsampling: interpolation method used  %
%                                              to assign a value to       %
%                                              every voxel of the 3D      %
%                                              volume after resampling    %
%           - ETL: echo train length, i.e. number of 180�-RF pulses       %
%           - flipAngle: refocusing flip angle (in degrees)               %
%           - ESP: echo spacing (in ms)                                   %
%           - TEeff: effective echo time (in ms)                          %
%           - TR: time from the application of an excitation pulse to     %
%                 the application of the next pulse (i.e., echo spacing   %
%                 in the EPG simulations) (in ms)                         %
%           - NbSlices: number of slices                                  %
%           - interleavedSlices_index: list of indexes of the slices      %
%                                      corresponding to an interleaved    %
%                                      acquisition scheme                 %
%           - Sl_to_Sl: slice profile, including slice gap if any         %
%           - SliceProfile: Gaussian slice profile                        %
%           - SliceThickness: thickness of a slice (in mm)                %
%           - FOVRead: dimension of the field-of-view in the read-out     %
%                      direction (in mm)                                  %
%           - FOVPhase: dimension of the field-of-view in the phase       %
%                       direction (in mm)                                 %
%           - BaseResolution: matrix size in the read-out direction (in   %
%                             voxels)                                     %
%           - nPE: number of phase encoding lines                         %
%           - ACF: acceleration factor                                    %
%           - RefLines: number of lines that are consecutively sampled    %
%                       around the center of K-space                      %
%           - GA: gestational age of the fetus (in weeks)                 %
%           - sub_id: subject ID                                          %
%           - WM_heterogeneity: boolean value that determines wether to   %
%                               implement (1) or not (0) white matter     %
%                               maturation processes in the simulated     %
%                               images                                    %
%           - affine: affine orientation matrix of the anatomical model   %
%                     of the fetal brain                                  %
%           - SimRes: resolution of the original 3D high-resolution       %
%                     anatomical model of the fetal brain (in mm)         %
%           - axcodes: axes codes of the 3D anatomical model of the       %
%                      fetal brain                                        %
%           - shift: displacement (in voxels) of the slice slab in the    %
%                    slice thickness direction                            %
%           - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
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
%  output:  - KSpace: Fourier domain of the simulated images              %
%           - Sl_Volume:          %
%                                                                         %
%                                                                         %
%  H�l�ne Lajous, 2023-03-17                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [   KSpace, ...
          Sl_Volume] = Kspace_sampling(                     T2decay_zp, ...
                                                           Fetal_Brain, ...
                                                           fetal_model, ...
                                                       b1map_upsampled, ...
                                                   Fetal_Brain_Tissues, ...
                                                             SimResReo, ...
                                                        SamplingFactor, ...
                                               motion_corrupted_slices, ...
                                              translation_displacement, ...
                                                        rotation_angle, ...
                                                         rotation_axis, ...
                                       interpolation_method_upsampling, ...
                                                                   ETL, ...
                                                             flipAngle, ...
                                                                   ESP, ...
                                                                 TEeff, ...
                                                                    TR, ...
                                                              NbSlices, ...
                                               interleavedSlices_index, ...
                                                              Sl_to_Sl, ...
                                                          SliceProfile, ...
                                                        SliceThickness, ...
                                                               FOVRead, ...
                                                              FOVPhase, ...
                                                        BaseResolution, ...
                                                                   nPE, ...
                                                                   ACF, ...
                                                              RefLines, ...
                                                                    GA, ...
                                                                sub_id, ...
                                                      WM_heterogeneity, ...
                                                                affine, ...
                                                                SimRes, ...
                                                           axcodes_reo, ...
                                                                 shift, ...
                                                           orientation, ...
                                                                 T1_WM, ...
                                                                 T2_WM, ...
                                                                 T1_GM, ...
                                                                 T2_GM, ...
                                                                T1_CSF, ...
                                                                T2_CSF)

% Input check
if nargin < 42
    error('Missing input(s).');
elseif nargin > 42
    error('Too many inputs.');
end

% Memory pre-allocation
Sl_Volume = zeros(size(T2decay_zp,1), size(T2decay_zp,2), NbSlices, size(T2decay_zp,4));

% Preload K-Space array which has dimensions of number of lines, number of
% Phase Encoding directions, and number of Slices
KSpace = zeros(BaseResolution, nPE, NbSlices, 'single');

% Create a "SLAB" array which is essentially the K-Space of our
% multislice volume
SLAB = zeros([BaseResolution, nPE, NbSlices, size(T2decay_zp,4)]);

SubunitRes = SimResReo(3) / SamplingFactor;  %mm

% Computation time to sample k-space
tic

motion_index = 0;

% Loop through slices
for iSlice=1:length(interleavedSlices_index)
    disp(['Slice ', num2str(interleavedSlices_index(iSlice)), ' of ', num2str(length(interleavedSlices_index))])
    index = interleavedSlices_index(iSlice)*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
    % Inter-slice motion -
    % - Random translation: uniform distribution between
    % [-translation_amplitude,+translation_amplitude]mm
    if any(interleavedSlices_index(iSlice)==motion_corrupted_slices)
        motion_index = motion_index + 1;
        Fetal_Brain_rotated = moving_3Dbrain(                             Fetal_Brain, ...
                                                                            SimResReo, ...
                                             translation_displacement(motion_index,:), ...
                                                                            'nearest', ...
                                                         rotation_angle(motion_index), ...
                                                        rotation_axis(motion_index,:), ...
                                                                            'nearest', ...
                                                                               'crop');
        Fetal_Brain_rotated_upsampled = sampling_OoP(            Fetal_Brain_rotated, ...
                                                                     SamplingFactor, ...
                                                     interpolation_method_upsampling);
        [ref_T1map_rotated, ref_T2map_rotated] = tissue_to_MR(                  fetal_model, ...
                                                              Fetal_Brain_rotated_upsampled, ...
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
                                                                                     T2_CSF);
        T2decay_moved = compute_t2decay(Fetal_Brain_rotated_upsampled, ...
                                                      b1map_upsampled, ...
                                                    ref_T1map_rotated, ...
                                                    ref_T2map_rotated, ...
                                                                  ETL, ...
                                                            flipAngle, ...
                                                                  ESP, ...
                                                      SamplingFactor);
        T2decay_moved_zp = Resize_Volume(T2decay_moved, size(T2decay_zp));
        T2decay_zp(:,:,index:index+round(SliceThickness/SubunitRes)-1,:) = T2decay_moved_zp(:,:,index:index+round(SliceThickness/SubunitRes)-1,:);
        % Update the T2decay_zp variable as it is unlikely that the
        % fetus comes back to its initial position at iSlice+1 and
        % following slices after he moved
        if iSlice < length(interleavedSlices_index)
            for nextSlice=iSlice+1:length(interleavedSlices_index)
                nextSlice_index = interleavedSlices_index(nextSlice)*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
                T2decay_zp(:,:,nextSlice_index:nextSlice_index+round(SliceThickness/SubunitRes)-1,:) = T2decay_moved_zp(:,:,nextSlice_index:nextSlice_index+round(SliceThickness/SubunitRes)-1,:);
            end
        end
        clear ref_T2map_rotated
        clear ref_T1map_rotated
        clear Fetal_Brain_rotated
        clear Fetal_Brain_rotated_upsampled
        clear T2decay_moved
        clear T2decay_moved_zp
    end
    % Multiply the 4D matrix by the slice profile
    T2decay_zp(:,:,index:index+round(SliceThickness/SubunitRes)-1,:) = T2decay_zp(:,:,index:index+round(SliceThickness/SubunitRes)-1,:).*permute(SliceProfile, [2,3,1]);
    % Sum the contribution of all voxels at the same (x,y) location across
    % the slice thickness direction
    Sl_Volume(:,:,interleavedSlices_index(iSlice),:) = sum(T2decay_zp(:,:,index:index+round(SliceThickness/SubunitRes)-1,:),3);
    % Simulation of k-space sampling as for FSE sequences
%     tic
    for iEcho=1:size(Sl_Volume,4)
        %clc;
%         disp(['Echo ', num2str(iEcho), ' of ', num2str(size(Sl_Volume,4))])
        SLAB(:,:,interleavedSlices_index(iSlice),iEcho) = Resize_Volume(fft2c(Resize_Volume(Sl_Volume(:,:,interleavedSlices_index(iSlice),iEcho), [round(FOVRead/SimResReo(1)), round(FOVPhase/SimResReo(2)), NbSlices])), [BaseResolution, nPE, NbSlices]);
    end
%     toc
%     tic
	% Calculating k-space sampling based on the desired echo time occuring
    % at the center of k-space
    if ACF~=1 && RefLines~=0
        temp = (nPE/2-RefLines/2)+RefLines+1:ACF:nPE;
        SamplingOrder = fliplr([ACF:ACF:(nPE/2-RefLines/2), (nPE/2-RefLines/2)+1:(nPE/2-RefLines/2)+RefLines, temp(1:round(TEeff/TR-RefLines/2))]);
        clear temp
    elseif ACF==1 && RefLines==0
        temp = nPE/2+1:ACF:nPE;
        SamplingOrder = fliplr([ACF:ACF:nPE/2, temp(1:round(TEeff/TR))]);
        clear temp
    elseif ACF~=1 && RefLines==0
        temp = nPE/2+ACF:ACF:nPE;
        SamplingOrder = fliplr([ACF:ACF:nPE/2, temp(1:round(TEeff/TR))]);
        clear temp
    end
    % Loop through phase encoding lines
    for iLine=1:length(SamplingOrder)
        KSpace(:, SamplingOrder(iLine), interleavedSlices_index(iSlice)) = SLAB(:, SamplingOrder(iLine), interleavedSlices_index(iSlice), iLine);
    end
    % This next block finds the missing lines and replaces them with the
    % closest echo time
    %Sum non-zero contributions in KSpace for slice 
    %interleavedSlices_index(iSlice), and echo iEcho
    if sum(sum(KSpace(:,:,interleavedSlices_index(iSlice))))~=0
        %The sampled lines in KSpace for slice iSlice, and echo iEcho
        %are non zero
        SampledLines = find(squeeze(KSpace(1,:,interleavedSlices_index(iSlice)))~=0);
        %Loop through all lines of KSpace
        for iLine=1:size(KSpace,2)
            %If a line was not sampled
            if KSpace(1,iLine,interleavedSlices_index(iSlice))==0
                %Find only the first non-zero line following iLine
                iFind = find(SampledLines>iLine,1,'first');
                %If a sampled line iFind was found following the
                %non-sampled line iLine, copy this iFind line from SLAB
                %to the corresponding line position in KSpace
                if ~isempty(iFind)
                    KSpace(:,iLine,interleavedSlices_index(iSlice)) = SLAB(:,iLine,interleavedSlices_index(iSlice),SamplingOrder==SampledLines(iFind));
                else
                    %If no sampled line was found following the non-sampled
                    %line iLine, use hermitian symmetry to fill KSpace: the
                    %line symmetrical to iLine compared to the center of
                    %KSpace is size(KSpace,2)-iLine+1
                    KSpace(:,iLine,interleavedSlices_index(iSlice)) = fliplr(conj(KSpace(:,size(KSpace,2)-iLine+1,interleavedSlices_index(iSlice)))')';
                end
            end
        end
    end
%     toc
end

% Display computation time
time3=toc;
fprintf('Computation time to sample k-space: %0.5f seconds.\n', time3);

end