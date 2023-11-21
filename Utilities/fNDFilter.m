%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   makeFilter      Creates an N-D filter                                 %
%                                                                         %
%     flt_kdat = fNDFilter(kdat, flt_name);                               %
%                                                                         %
%     INPUT:    kdat        N-D data to filter                            %
%               flt_name    Filter name to use                            %
%               dims        Array of dimensions to filter                 %
%                                                                         %
%     OUTPUT:   flt_kdat    Filtered N-D data                             %
%                                                                         %
%                                                                         %
%  (c) Jerome Yerly, 2013                                                 %
%  yerlyj.mri@gmail.com                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function flt_kdat = fNDFilter(kdat, flt_name, dims)
  
  % Check input and output arguments
    narginchk(2,3); % Between two and three input arguments
    nargoutchk(1,1); % 1 output agrument

    ndim     = ndims(kdat); % Number of dimensions
    
    if nargin < 3
        dims = 1:ndim;   % Filter all dimensions by default
    else
      % Verify that the array dims contains feasible values, i.e. between 1 and ndim
        if min(dims) < 1 && max(dims) > ndim
            error('The array dims can only contain values between 1 and %d. min(dims)= %d and max(dims)= %d',ndim,max(dims),min(dims));
        end
    end

    dims     = unique(dims); % Make sure that values in dims are unique
    nfilter  = length(dims); % Number of dimensions to filter
    
  % Create a distance matrix
    x = cell(ndim,1);
    for i = 1:ndim
        x{i} = linspace( -1, 1, size(kdat,i) );
    end
    X = meshgrid_ND(x);
  
    kr = zeros(size(X{1}));
    for i = 1:nfilter
        dim = dims(i);
        kr = kr + X{dim}.^2;
    end
    kr = sqrt(kr);
        
  % Define filter
    if ( strcmpi(flt_name,'von Hann') == 1 ),
      tmp_alf = 0.50;
      flt = (1 - tmp_alf) + tmp_alf * cos( pi * kr );
    elseif ( strcmpi(flt_name,'Hamming') == 1 ),
      tmp_alf = 0.46;
      flt = (1 - tmp_alf) + tmp_alf * cos( pi * kr );
    elseif ( strcmpi(flt_name,'Fermi') == 1 ),
      % Suggested values:
      %   - tmp_fR = 0.85
      %   - tmp_fW = 23
      tmp_fR  = 0.85;   % cutoff at 50% max
      tmp_fW  = 23;     % 1/transition width
      flt = 1 ./ ( 1 + exp( -tmp_fW * (tmp_fR - kr) ) );
    else
      flt = ones ( size ( kdat ) );
    end
    
  % Filter k-space
    flt_kdat = kdat .* flt;

  % Return
    return

end