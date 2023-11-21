%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that reads the segmented high-resolution anatomical MR        %
%  images of the fetal brain from which the simulated images will be      %
%  derived at a given gestational age (GA).                               %
%  Here, we consider segmented high-resolution images from Gholipour, A.  %
%  et al. A normative spatiotemporal MRI atlas of the fetal brain for     %
%  automatic segmentation and analysis of early brain growth. Scientific  %
%  Reports 7, 476 (2017). https://doi.org/10.1038/s41598-017-00525-w      %
%                                                                         %
%                ornt = io_orientation(affine, tolerance);                %
%                                                                         %
%  inputs:  - affine: transformation affine matrix from the anatomical    %
%                     fetal brain model                                   %
%           - tolerance: threshold below which SVD values of the affine   %
%                        are considered zero                              %
%                                                                         %
%  output:  - ornt: array of the closest output axis (R,A,S) and          %
%                   direction (1 if the input axis is in the same         %
%                   direction as the corresponding output axis and -1 if  %
%                   it is in the opposite direction                       %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2023-01-30                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ornt = io_orientation(affine, varargin)

% Input check
if nargin < 1
    error('Missing input(s).');
elseif nargin==1
    tolerance = "None"; %default value for the tolerance threshold
else
    if numel(varargin) > 0  %optional input arguments are provided
        if numel(varargin) < 2
            error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
        end
        switch varargin{1}
            case 'tolerance'
                tolerance = varargin{2};
            otherwise
                error('Unexpected ''ParameterName'' input: %s\n', varargin{1});
        end
    end
end

% Extract the underlying rotation, zoom, shear matrix
RZS = affine(1:size(affine,1)-1, 1:size(affine,2)-1);
zooms = sqrt(sum(RZS .* RZS, 1));
% Zooms can be zero, in which case all elements in the column are zero, and
% we can leave them as they are
zooms(zooms==0) = 1;
RS = RZS ./ zooms;
% Transform below is polar decomposition, returning the closest shearless
% matrix R to RS
[P, S, Qs] = svd(RS);
% Threshold the singular values to determine the rank
if tolerance=="None"
    tolerance = max(S) * max(size(RS)) * eps(class(S));
end
keep = svd(RS)' > tolerance;
R = P(:,keep) * Qs(keep,keep);

% The matrix R is such that the scalar product between R and R*T is 
% projection onto the columns of P(:,keep) and the scalar product between 
% R.T and R is the projection onto the rows of Qs(keep). R gives rotation 
% of the unit input vectors to output coordinates. Therefore, the row index
% of abs max R(:,N) is the output axis changing most as input axis N 
% changes. In case there are ties, we choose the axes iteratively, removing
% used axes from consideration as we go.
ornt = NaN(size(affine,2)-1,2);
for in_ax=1:size(affine,2)-1
    col = R(:,in_ax);
    if any(col)
        [Max,indexMax] = max(abs(col));
        out_ax = indexMax;
        ornt(in_ax,1) = out_ax;
        assert(col(out_ax)~=0)
        if col(out_ax)<0
            ornt(in_ax,2) = -1;
        else
            ornt(in_ax,2) = 1;
        end
        % Remove the identified axis from further consideration, by zeroing
        % out the corresponding row in R
        R(out_ax, :) = 0;
    end
end

end