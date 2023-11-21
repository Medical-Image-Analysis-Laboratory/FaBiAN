%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that resamples a 3D volume by a given sampling factor and     %
%  according to the specified interpolation method in the slice           %
%  thickness direction (3. index, "Out-of-Plane").                        %
%                                                                         %
%           VolumeResampled = sampling_OoP(             Volume, ...       %
%                                               SamplingFactor, ...       %
%                                          InterpolationMethod);          %
%                                                                         %
%  inputs:  - Volume: 3D volume to be resampled                           %
%           - SamplingFactor: factor of up- or down-sampling              %
%           - InterpolationMethod: interpolation method used to assign a  %
%                                  value to every voxel of the resampled  %
%                                  3D volume. Images generated for the    %
%                                  qualitative evaluation by the          %
%                                  radiologists were simulated from       %
%                                  partial volume maps bilinearly         %
%                                  interpolated. To reduce the            %
%                                  computational burden that arises from  %
%                                  EPG simulations with many non-unique   %
%                                  combinations of (b1,T1,T2), the        %
%                                  images generated for the data          %
%                                  augmentation experiment were           %
%                                  simulated from partial volume maps     %
%                                  interpolated using a nearest-          %
%                                  -neighboor method.                     %
%                                                                         %
%  output:  - VolumeResampled: 3D volume after resampling in the slice    %
%                               thickness direction                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function VolumeResampled = sampling_OoP(             Volume, ...
                                             SamplingFactor, ...
                                        InterpolationMethod)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs.');
end


switch InterpolationMethod
    case 'nearest'
        for index=1:size(Volume,3)
            upsampling = index*SamplingFactor-(SamplingFactor-1);
            VolumeResampled(:,:,upsampling:upsampling+SamplingFactor-1) = repmat(Volume(:,:,index), [1, 1, SamplingFactor]);
        end
    case 'linear'
        Volume_orient = reshape(Volume, [size(Volume,1)*size(Volume,2), size(Volume,3)]);
        Volume_up = imresize(Volume_orient, [size(Volume_orient,1) size(Volume_orient,2)*SamplingFactor], 'bilinear');
        VolumeResampled = reshape(Volume_up, [size(Volume,1), size(Volume,2), size(Volume_up,2)]);
end

% Display message for debugging
sprintf('Volume resampled.')

end