%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fZeroPadnDArray      Zero-pads an N-D array                           %
%                                                                         %
%     arrayOut = fZeroPadnDArray(arrayIn, newN);                          %
%                                                                         %
%                                                                         %
%  (c) Jerome Yerly 2013                                                  %
%  yerlyj.mri@gmail.com                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function arrayOut = fZeroPadnDArray( arrayIn, newN )

  % Number of dimension in arrayIn
    Ndim = ndims(arrayIn);
    
  % newN should be an array of length equal to the number of dimension in
  % arrayIn
    if length(newN) ~= Ndim
        return;
    end

  % Get dimension of array to pad
    oldN = size( arrayIn );
    
    if isreal(arrayIn)
      % Initialize the padded array
        arrayOut = zeros( newN, class(arrayIn) );
    else
      % Initialize the padded array
        arrayOut = complex(zeros( newN, class(arrayIn) ));
    end
    
  % Calculate the old and new indice corresponding to the centre of the
  % data
    idxCtr_old = floor( oldN / 2 ) + 1;
    idxCtr_new = floor( newN / 2 ) + 1;

  % Calculate the start indices where to place the centered data
    idxStart = idxCtr_new - idxCtr_old + 1;
    
  % Calculate the end indices where to place the centered data
    idxEnd = idxStart + oldN - 1;
    
  % Loop through all dimension to generate indices
    for i = 1:Ndim
        
      % Construct command
        str = sprintf('idx%d = idxStart(%d) : idxEnd(%d);', i, i, i);
        
      % Eval command
        eval(str);
        
    end
    
  % Initialize string for command
    tmpStr = '';
    
  % Loop through all dimension to construct the next command
    for i = 1:Ndim
        
        if i == 1
          % Construct string for next command
            tmpStr = sprintf('idx%d ', i);
        else
          % Construct string for next command
            tmpStr = sprintf('%s, idx%d ', tmpStr, i);
        end
    end
    
  % Construct command
    str = sprintf('arrayOut(%s) = arrayIn;', tmpStr);
    
  % Eval command
    eval(str);

end