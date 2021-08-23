%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that reads the high-resolution segmented anatomical MR        %
%  images of the fetal brain from which the simulated images will be      %
%  derived.                                                               %
%  Here, we consider high-resolution segmented images from Gholipour, A.  %
%  et al. A normative spatiotemporal MRI atlas of the fetal brain for     %
%  automatic segmentation and analysis of early brain growth. Scientific  %
%  Reports 7, 1–13, DOI: 10.1038/s41598-017-00525-w (2017).               %
%                                                                         %
%                   Fetal_Brain = brain_model(path, GA)                   %
%                                                                         %
%  input:   - path: folder where the high-resolution segmented images of  %
%                   the fetal brain are stored                            %
%           - GA: gestational age of the fetus (in weeks)                 %
%                                                                         %
%  output:  - Fetal_Brain: 3D volume of the fetal brain                   %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-30                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Fetal_Brain = brain_model(path, ...
                                     GA)


% Input check
if nargin < 2
    error('Missing input(s).');
elseif nargin > 2
    error('Too many inputs.');
end


% Load segmented HR images of the fetal brain at gestational age GA
Fetal_Brain = niftiread(strcat(path, 'STA', sprintf('%02s', num2str(GA)), '_tissue.nii.gz'));
