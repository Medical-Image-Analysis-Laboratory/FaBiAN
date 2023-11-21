%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that returns the axis direction codes for a given affine      %
%  matrix.                                                                %
%                                                                         %
%            axcodes = aff2axcodes(affine, labels, tolerance);            %
%                                                                         %
%  inputs:  - affine: affine transformation matrix from the anatomical    %
%                     fetal brain model                                   %
%           - labels: labels for negative and positive ends of output     %
%                     axes of the affine matrix                           %
%           - tolerance: threshold below which SVD values of the affine   %
%                        are considered zero                              %
%                                                                         %
%  output:  - axcodes: labels for positive end of voxel axes. Dropped     %
%                      axes get a label of None.                          %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2023-01-30                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function axcodes = aff2axcodes(affine, varargin)

% Input check
if nargin < 1
    error('Missing input(s).');
else
    defaultLabels = {["L","R"], ["P","A"], ["I","S"]};
    defaultTolerance = "None";
    p = inputParser;
    addRequired(p, 'affine');
    addOptional(p, 'labels', defaultLabels);
    addOptional(p, 'tolerance', defaultTolerance);
    parse(p, affine, varargin{:});
    
    ornt = io_orientation(p.Results.affine, 'tolerance', p.Results.tolerance);
    axcodes = ornt2axcodes(ornt, 'labels', p.Results.labels);
end

end