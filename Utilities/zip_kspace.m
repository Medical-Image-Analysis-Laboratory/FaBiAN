%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that filters data in K-space to avoid Gibbs ringing           %
%  artifacts and that performs zero-interpolation filling (ZIP), namely   %
%  fills the edges of the acquired K-space with zeros to reach the        %
%  desired reconstruction matrix size.                                    %
%                                                                         %
%              KSpace_zip = zip_kspace(KSpace, reconMatrix);              %
%                                                                         %
%  inputs   - KSpace: Fourier domain of the simulated images              %
%           - reconMatrix: size of the reconstruction matrix (in voxels)  %
%                                                                         %
%  output:  - KSpace_zip: Fourier domain of the simulated images after    %
%                         ZIP                                             %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-05-05                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function KSpace_zip = zip_kspace(     KSpace, ...
                                 reconMatrix)

% Input check
if nargin < 2
    error('Missing input(s).');
elseif nargin > 2
    error('Too many inputs.');
end

% Filter data in k-space to avoid Gibbs ringing artifacts
KSpace_filt = fNDFilter(KSpace, 'Fermi', [1,2]);

% Interpolate data via zero padding in k-space
KSpace_zip = fZeroPadnDArray(KSpace_filt, reconMatrix);

end