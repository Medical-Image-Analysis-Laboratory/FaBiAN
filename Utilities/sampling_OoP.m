%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that resamples a 3D volume by a given sampling factor and     %
%  according to the specified interpolation method in the slice           %
%  thickness direction (3. index).                                        %
%                                                                         %
%          Volume_resampled = sampling_OoP(              Volume, ...      %
%                                               sampling_factor, ...      %
%                                          interpolation_method)          %
%                                                                         %
%  input:   - Volume: 3D volume to be resampled                           %
%           - sampling_factor: factor of (up-)/(down-)sampling            %
%           - interpolation_method: method to be used for the             %
%                                   interpolation                         %
%                                                                         %
%  output:  - Volume_resampled: resampled 3D volume                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
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

% Display message for debugging
sprintf('Volume resampled.')
