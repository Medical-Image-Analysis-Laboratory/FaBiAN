%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that resamples a 3D volume in-plane                           %
%                                                                         %
%  inputs:  - Volume: 3D volume to be resampled                           %
%           - resampling_read: resampling in the readout direction        %
%           - resampling_phase: resampling in the phase-encoding          %
%                               direction                                 %
%           - interpolation_method: interpolation method used to          %
%                                   resample the volume, linear for an    %
%                                   image, nearest neighbors for a label  %
%                                   map or a binary mask                  %
%                                                                         %
%  output:  - Resampled_Volume: 3D volume after resampling in the         %
%                               in-plane orientation                      %
%                                                                         %
%  Hélène Lajous, 2021-07-05                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ResampledVolume = resampling_inplane(             Volume, ...
                                                   ResamplingRead, ...
                                                  ResamplingPhase, ...
                                              InterpolationMethod)


% Input check
if nargin < 4
    error('Missing input(s).');
elseif nargin > 4
    error('Too many inputs.');
end


switch InterpolationMethod
    case 'nearest'
        ResampledVolume = imresize(Volume, [round(size(Volume,1)/ResamplingRead) round(size(Volume,2)/ResamplingPhase)], 'nearest');
    case 'linear'
        ResampledVolume = imresize(Volume, [round(size(Volume,1)/ResamplingRead) round(size(Volume,2)/ResamplingPhase)], 'bilinear');
end

% Display message for debugging
sprintf('Volume resampled.')
