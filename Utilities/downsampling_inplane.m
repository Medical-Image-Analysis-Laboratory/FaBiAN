%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that downsamples a 3D volume in-plane according to the        %
%  interpolation method given:                                            %
%     - linear if it's an image                                           %
%     - nearest neighbors if it's a label map or a binary mask            %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-07-05                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Volume_downsampled = downsampling_inplane(              Volume, ...
                                                      downsampling_read, ...
                                                     downsampling_phase, ...
                                                   interpolation_method)


% Input check
if nargin < 4
    error('Missing input(s).');
elseif nargin > 4
    error('Too many inputs.');
end


switch interpolation_method
    case 'nearest'
        Volume_downsampled = imresize(Volume, [size(Volume,1)/downsampling_read size(Volume,2)/downsampling_phase], 'nearest');
    case 'linear'
        Volume_downsampled = imresize(Volume, [size(Volume,1)/downsampling_read size(Volume,2)/downsampling_phase], 'bilinear');
end

% Display message for debugging
sprintf('Volume downsampled.')
