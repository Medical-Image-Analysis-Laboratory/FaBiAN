%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   meshgrid_ND      Creates rectangular grid in N-D space                %
%                                                                         %
%     X = meshgrid_ND(x);                                                 %
%     Replicates the grid vectors in cell array x to produce a full grid  %
%                                                                         %
%     INPUT:    x     Cell array containing the grid vectors (N vectors)  %
%                                                                         %
%     OUTPUT:   X     Cell array containing the output coordinate arrays  %
%                     of dimension N                                      %
%                                                                         %
%                                                                         %
%  (c) Jerome Yerly 2013                                                  %
%  yerlyj.mri@gmail.com                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function X = meshgrid_ND(x)

ndim     = length(x); % Number of dimension
X        = cell(ndim,1);% Create cell array to contain all the grids
dimSizes = zeros(ndim,1); % Create vector to contain thd dimension of each dimension

% Get size of each dimension
for i = 1:ndim
    dimSizes(i) = numel(x{i});
end

% Create grids
for i = 1:ndim
    
    % Make sure x is a vector along dimension i.
    xx          = x{i};
    vecShape    = ones(1,ndim);
    vecShape(i) = dimSizes(i);
    xx          = reshape(xx(:),vecShape);
    
    % Repeat vector to create grid
    tmp     = dimSizes;
    tmp(i)  = 1;
    
    if ~isscalar(xx)
        xx      = repmat(xx,tmp(:)');
    end
    
    % Store in cell array
    X{i}    = xx;
end

end