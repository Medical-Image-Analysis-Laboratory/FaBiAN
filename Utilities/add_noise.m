%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that adds complex Gaussian noise (mean=0, standard deviation  %
%  to be defined by the user) to a volume in the Fourier domain.          %
%                                                                         %
%               KSpace_noise = add_noise(KSpace, std_noise);              %
%                                                                         %
%  inputs:  - KSpace: Fourier domain of the simulated images              %
%           - std_noise: standard deviation of the Gaussian noise added   %
%                        to the Fourier domain of the simulated images    %
%                                                                         %
%  output:  - KSpace_noise: Fourier domain of the simulated images after  %
%                           addition of complex Gaussian noise            %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-22                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function KSpace_noise = add_noise(KSpace, std_noise)

% Input check
if nargin < 2
    error('Missing input(s).');
elseif nargin > 2
    error('Too many inputs.');
end

cgnoise = complex(std_noise.*randn(size(KSpace)), std_noise.*randn(size(KSpace)));
KSpace_noise = KSpace + cgnoise;

end