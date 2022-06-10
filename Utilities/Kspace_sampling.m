%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that samples the Fourier domain (K-space) of the simulated    %
%  images according to the acquisition scheme implemented in Fast Spin    %
%  Echo (FSE) sequences.                                                  %
%                                                                         %
%        KSpace = Kspace_sampling(                     T2decay_zp, ...    %
%                                                 b1map_upsampled, ...    %
%                                             Fetal_Brain_Tissues, ...    %
%                                                          SimRes, ...    %
%                                                 sampling_factor, ...    %
%                                                       
%                                                              B0, ...    %
%                                         motion_corrupted_slices, ...    %
%                                        translation_displacement, ...    %
%                                                  rotation_angle, ...    %
%                                                   rotation_axis, ...    %
%                                 interpolation_method_upsampling, ...    %
%                                                             ETL, ...    %
%                                                             ESP, ...    %
%                                                           TEeff, ...    %
%                                                              TR, ...    %
%                                                        NbSlices, ...    %
%                                         interleavedSlices_index, ...    %
%                                                        Sl_to_Sl, ...    %
%                                                    SliceProfile, ...    %
%                                                  SliceThickness, ...    %
%                                                         FOVRead, ...    %
%                                                        FOVPhase, ...    %
%                                                  BaseResolution, ...    %
%                                                             nPE, ...    %
%                                                             ACF, ...    %
%                                                       RefLines);        %
%                                                                         %
%  inputs:  - T2decay_zp: 4D matrix that combines the anatomical          %
%                         information from the segmented high-resolution  %
%                         fetal brain volume and the T2 decay computed    %
%                         in every voxel of this volume after zero-       %
%                         -padding                                        %
%           - Fetal_Brain: segmented high-resolution 3D volume of the     %
%                          fetal brain                                    %
%           - b1map_upsampled: B1+ bias field map resampled by the        %
%                              sampling_factor in the slice thickness     %
%                              direction                                  %
%           - Fetal_Brain_Tissues: list of tissues in the fetal brain     %
%                                  that were segmented and labeled        %
%           - SimRes: resolution of the original high-resolution fetal    %
%                     brain images (isotropic, in mm)                     %
%           - sampling_factor: factor of up- or down-sampling             %
%           - B0: main magnetic field strength                            %
%           - motion_corrupted_slices: list of indexes of the slices      %
%                                      corrupted by motion                %
%           - translation_displacement: table of translation              %
%                                       displacements (in mm) along the   %
%                                       three main axes for each motion-  %
%                                       -corrupted slice                  %
%           - rotation_angle: list of rotation angles (in °) for 3D       %
%                             rotation in every motion-corrupted slice    %
%           - rotation_axis: table of rotation axes for 3D rotation in    %
%                            every motion-corrupted slice                 %
%           - interpolation_method_upsampling: interpolation method used  %
%                                              to assign a value to       %
%                                              every voxel of the 3D      %
%                                              volume after resampling    %
%           - ETL: echo train length, i.e. number of 180°-RF pulses       %
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
%                                                                         %
%  output:  - KSpace: Fourier domain of the simulated images              %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function KSpace = Kspace_sampling(                     T2decay_slice, ...    
                                                           SimRes, ...
                                                  sampling_factor, ...
                                                            TEeff, ...
                                                               TR, ...
                                                            index, ...
                                                     SliceProfile, ...
                                                   SliceThickness, ...
                                                          FOVRead, ...
                                                         FOVPhase, ...
                                                   BaseResolution, ...
                                                              nPE, ...
                                                              ACF, ...
                                                         RefLines)

% Input check
if nargin < 14
    error('Missing input(s).');
elseif nargin > 14
    error('Too many inputs.');
end

% Memory pre-allocation
Sl_Volume = zeros(size(T2decay_slice,1), size(T2decay_slice,2), 1, size(T2decay_slice,4));

% Preload K-Space array which has dimensions of number of lines, number of
% Phase Encoding directions, and number of Slices
KSpace = zeros(BaseResolution, nPE, 1, 'single');

% Create a "SLAB" array which is essentially the K-Space of our
% multislice volume
SLAB = zeros([BaseResolution, nPE, 1, size(T2decay_slice,4)]);

SubunitRes = SimRes / sampling_factor;  %mm



% Multiply the 4D matrix by the slice profile
T2decay_slice = T2decay_slice.*permute(SliceProfile, [2,3,1]);
% Sum the contribution of all voxels at the same (x,y) location across
% the slice thickness direction
Sl_Volume = sum(T2decay_slice,3);
% Simulation of k-space sampling as for FSE sequences

for iEcho=1:size(Sl_Volume,4)
    clc;
    disp(['Echo ', num2str(iEcho), ' of ', num2str(size(Sl_Volume,4))])
    SLAB(:,:,: ,iEcho) = Resize_Volume(fft2c(Resize_Volume(Sl_Volume(:,:,:,iEcho), [FOVRead/SimRes, round(FOVPhase)/SimRes, 1])), [BaseResolution, nPE, 1]);
end

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
end
% Loop through phase encoding lines
for iLine=1:length(SamplingOrder)
    KSpace(:, SamplingOrder(iLine),:) = SLAB(:, SamplingOrder(iLine),:, iLine);
end
% This next block finds the missing lines and replaces them with the
% closest echo time
%Sum non-zero contributions in KSpace for slice 
%interleavedSlices_index(iSlice), and echo iEcho
if sum(sum(KSpace))~=0
    %The sampled lines in KSpace for slice iSlice, and echo iEcho
    %are non zero
    SampledLines = find(squeeze(KSpace(1,:,:))~=0);
    %Loop through all lines of KSpace
    for iLine=1:size(KSpace,2)
        %If a line was not sampled
        if KSpace(1,iLine,:)==0
            %Find only the first non-zero line following iLine
            iFind = find(SampledLines>iLine,1,'first');
            %If a sampled line iFind was found following the
            %non-sampled line iLine, copy this iFind line from SLAB
            %to the corresponding line position in KSpace
            if ~isempty(iFind)
                KSpace(:,iLine,:) = SLAB(:,iLine,:,SamplingOrder==SampledLines(iFind));
            else
                %If no sampled line was found following the non-sampled
                %line iLine, use hermitian symmetry to fill KSpace: the
                %line symmetrical to iLine compared to the center of
                %KSpace is size(KSpace,2)-iLine+1
                KSpace(:,iLine,:) = fliplr(conj(KSpace(:,size(KSpace,2)-iLine+1,:))')';
            end
        end
    end
end


% Display computation time
fprintf('Computation time to sample k-space: %0.5f seconds.\n', toc);

end