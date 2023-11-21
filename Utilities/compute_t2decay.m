%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that computes the T2 decay in every voxel of the fetal brain  %
%  volume from reference T1 and T2 maps, sequence parameters and          %
%  intensity non-uniformity fields, following the Extended Phase Graph    %
%  (EPG) formalism.                                                       %
%  The EPG concept is described in details in: Weigel, M. Extended phase  %
%  graphs: Dephasing, RF pulses, and echoes - pure and simple. Journal    %
%  of Magnetic Resonance Imaging 41, 266-295 (2015).                      %
%  https://doi.org/10.1002/jmri.24619. The associated EPG simulation      %
%  code from Matthias Weigel for multi-spin echo sequences                %
%  (cp_cpmg_epg_domain_fplus_fminus.m) used in the following can be       %
%  downloaded here: https://github.com/matthias-weigel/EPG.               %
%                                                                         %
%            T2decay = compute_t2decay(Fetal_Brain_upsampled, ...         %
%                                            b1map_upsampled, ...         %
%                                                  ref_T1map, ...         %
%                                                  ref_T2map, ...         %
%                                                        ETL, ...         %
%                                                  flipAngle, ...         %
%                                                        ESP, ...         %
%                                            sampling_factor);            %
%                                                                         %
%  inputs:  - Fetal_Brain_upsampled: segmented high-resolution 3D volume  %
%                                    of the fetal brain after upsampling  %
%                                    by a given sampling factor in the    %
%                                    slice thickness direction            %
%           - b1map_upsampled: 3D B1+ bias field map after upsampling by  %
%                              a given sampling factor in the slice       %
%                              thickness direction                        %
%           - ref_T1map: reference T1 map of the fetal brain (3D)         %
%           - ref_T2map: reference T2 map of the fetal brain (3D)         %
%           - ETL: echo train length, i.e. number of 180°-RF pulses       %
%           - flipAngle: refocusing flip angle (in degrees)               %
%           - ESP: echo spacing (in ms)                                   %
%           - sampling_factor: factor by which the fetal brain volume     %
%                              and the B1 bias field have been upsampled  %
%                                                                         %
%  output:  - T2decay: 4D matrix that combines the anatomical             %
%                      information from the segmented high-resolution     %
%                      fetal brain volume and the T2 decay computed in    %
%                      every voxel of this volume                         %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-20                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function T2decay = compute_t2decay(Fetal_Brain_upsampled, ...
                                         b1map_upsampled, ...
                                               ref_T1map, ...
                                               ref_T2map, ...
                                                     ETL, ...
                                               flipAngle, ...
                                                     ESP, ...
                                         sampling_factor)

% Input check
if nargin < 8
    error('Missing input(s).');
elseif nargin > 8
    error('Too many inputs.');
end

% Unwrap maps for computation speed
unwrap_b1map = b1map_upsampled(:);
unwrap_ref_T1map = ref_T1map(:);
unwrap_ref_T2map = ref_T2map(:);

% Find every possible combination of {b1, T1, T2} values (umap) and store
% the corresponding row where this combination first occurs (imap)
[umap, imap] = unique([unwrap_b1map, unwrap_ref_T1map, unwrap_ref_T2map], 'rows');
% Exclude background
imap = imap(umap(:,2)~=0, :);
umap = umap(umap(:,2)~=0, :);

uT2decay = zeros(size(umap,1), ETL);

% Computation time
tic

parfor i=1:size(uT2decay,1)
    uT2decay(i, :) = real(cp_cpmg_epg_domain_fplus_fminus(umap(i,1).*90, ETL, umap(i,1).*flipAngle, ESP, umap(i,2), umap(i,3)))/sampling_factor;
end

T2decay = zeros(length(Fetal_Brain_upsampled(:)), ETL, 'single');
for i=1:size(uT2decay,1)
    j = find(unwrap_b1map==umap(i,1) & unwrap_ref_T1map==umap(i,2) & unwrap_ref_T2map==umap(i,3));
    for k=1:length(j)
        T2decay(j(k),:) = uT2decay(i,:);
    end
end
T2decay = reshape(T2decay, [size(Fetal_Brain_upsampled), size(T2decay,2)]);

% Display computation time
time2=toc;
fprintf('Computation time to run EPG simulations in every voxel of the image: %0.5f seconds.\n', time2);

end