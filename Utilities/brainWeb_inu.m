%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that calls simulated intensity non-uniformity fields (3D),    %
%  resizes them to the fetal brain volume dimensions and normalizes them  %
%  to a level of 40% (i.e., to a range of values of 0.8-1.2 over the      %
%  fetal brain volume).                                                   %
%                                                                         %
%                  b1map = brainWeb_inu(inu, Fetal_Brain)                 %
%                                                                         %
%  input:   - inu: simulated intensity non-uniformity fields (3D). Here,  %
%                  we consider B1 bias field variations in raw byte       %
%                  (unsigned) format obtained from the BrainWeb database  %
%                  -- BrainWeb: Simulated brain database.                 %
%                  https://brainweb.bic.mni.mcgill.ca/brainweb/ --        %
%           - Fetal_Brain: fetal brain volume                             %
%                                                                         %
%  output:  - b1map: simulated B1 bias field variations. This 3D volume   %
%                    is the same size as the fetal brain volume.          %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2020-09-30                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function b1map = brainWeb_inu(inu, Fetal_Brain)

% Input check
if nargin < 2
    error('Missing input(s).');
elseif nargin > 2
    error('Too many inputs.');
end

% x, y and z are the dimensions of the simulated intensity non-uniformity
% fields from the BrainWeb database
z = 181;
y = 217;
x = 181;
volume = zeros(x,y,z);

% Read BrainWeb intensity non-uniformity fields (3D)
f = fopen(inu, 'r');
for i=1:z
    Im = fread(f, [x,y], 'uint8');
    volume(:,:,i) = Im;
end

% Resize the intensity non-uniformity fields to match the fetal brain
% volume dimensions
b1map_res = Resize_Volume(volume, size(Fetal_Brain));

% Normalize the intensity non-uniformity fields by 1.2 (i.e., level of 40%)
b1map = b1map_res ./ max(max(max(b1map_res))) * 1.2;

fclose(f);
