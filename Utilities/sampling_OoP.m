%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that resamples a 3D volume by a given sampling factor and     %
%  according to the specified interpolation method in the slice           %
%  thickness direction (3. index, "Out-of-Plane").                        %
%                                                                         %
%          Volume_resampled = sampling_OoP(              Volume, ...      %
%                                               sampling_factor, ...      %
%                                          interpolation_method);         %
%                                                                         %
%  inputs:  - Volume: 3D volume to be resampled                           %
%           - sampling_factor: factor of up- or down-sampling             %
%           - interpolation_method: interpolation method used to assign   %
%                                   a value to every voxel of the 3D      %
%                                   volume after resampling               %
%                                                                         %
%  output:  - Volume_resampled: 3D volume after resampling in the slice   %
%                               thickness direction                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Volume_resampled = sampling_OoP(              Volume, ...
                                              sampling_factor, ...
                                         interpolation_method)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs.');
end

if length(sampling_factor)==1

    switch interpolation_method
        case 'nearest'
            for index=1:size(Volume,3)
                upsampling = index*sampling_factor-(sampling_factor-1);
                Volume_resampled(:,:,upsampling:upsampling+sampling_factor-1) = repmat(Volume(:,:,index), [1, 1, sampling_factor]);
            end
        case 'linear'
            Volume_orient = reshape(Volume, [size(Volume,1)*size(Volume,2), size(Volume,3)]);
            Volume_up = imresize(Volume_orient, [size(Volume_orient,1) size(Volume_orient,2)*sampling_factor], 'bilinear');
            Volume_resampled = reshape(Volume_up, [size(Volume,1), size(Volume,2), size(Volume_up,2)]);
    end
elseif length(sampling_factor)==3
    Volume_resampled = imresize3(Volume, size(Volume).*sampling_factor, 'bilinear');
    
end
% Display message for debugging
sprintf('Volume resampled.')

end