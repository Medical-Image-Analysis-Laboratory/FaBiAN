%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that reorients a 3D volume to have its orientation, which     %
%  corresponds to the slice thickness direction, given by the 3. index.   %
%                                                                         %
%        Volume_Reoriented = volume_reorient(Volume, orientation);        %
%                                                                         %
%  inputs:  - Volume: 3D volume to be reoriented                          %
%           - Affine: affine transformation matrix from the anatomical    %
%                     fetal brain model                                   %
%           - SimRes: resolution of the original 3D anatomical model of   %
%                     the fetal brain                                     %
%           - out_axcodes: desired output axes codes                      %
%                                                                         %
%  output:  - ReorientedVolume: 3D volume after reorientation, i.e. its   %
%                               third dimension corresponds to the slice  %
%                               thickness direction                       %
%           - NewAffine: affine transformation matrix of the reoriented   %
%                        fetal brain model                                %
%           - NewSimRes: resolution of the reoriented 3D anatomical       %
%                        model of the fetal brain                         %
%                                                                         %
%                                                                         %
%  Oscar Esteban, Hélène Lajous, 2023-02-27                               %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ReorientedVolume, ...
                 NewAffine, ...
                 NewSimRes] = reorient_volume(     Volume, ...
                                                   Affine, ...
                                                   SimRes, ...
                                              out_axcodes)

% Input check
if nargin < 4
    error('Missing input(s).');
elseif nargin > 4
    error('Too many inputs.');
end

% Axes codes of the input volume
in_axcodes = aff2axcodes(Affine);

% Create a dictionary with the axes codes of the original model
codemap = containers.Map({char(in_axcodes(1)) char(in_axcodes(2)) char(in_axcodes(3))}, [1 2 3]);

for i=1:length(out_axcodes)
    ornt(i) = codemap(out_axcodes(i));
end

% Reorient volume according to the desired output axes codes and update the
% corresponding affine
ReorientedVolume = permute(Volume, ornt);
NewAffine = Affine(:, [ornt 4]);
% NewRotations = Affine(1:3,ornt)~=0;
% NewRotations(4,4) = 1;
% NewAffine = NewRotations * Affine;
% NewAffine = double(Affine(1:3,ornt)~=0);
% NewAffine(4,4) = 1;
% NewRotations = 
% NewAffine(:,4) = NewAffine * Affine(:,4);
NewSimRes = SimRes(ornt);

% Display message for debugging
sprintf('The volume was successfully reoriented.')

end