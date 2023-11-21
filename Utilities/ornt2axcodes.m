%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that converts orientation to labels for axis directions.      %
%                                                                         %
%                  axcodes = ornt2axcodes(ornt, labels);                  %
%                                                                         %
%  inputs:  - ornt: orientation array corresponding to the anatomical     %
%                   model of the fetal brain                              %
%           - labels: (begininng, end) of output axis                     %
%                                                                         %
%  output:  - axcodes: labels for positive end of voxel axes.             %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2023-01-30                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function axcodes = ornt2axcodes(ornt, varargin)

% Input check
if nargin < 1
    error('Missing input(s).');
elseif nargin==1
    labels = {["L","R"], ["P","A"], ["I","S"]};   %default value for labels
else
    if numel(varargin) > 0  %optional input arguments are provided
        if numel(varargin) < 2
            error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
        end
        switch varargin{1}
            case 'labels'
                labels = varargin{2};
            otherwise
                error('Unexpected ''ParameterName'' input: %s\n', varargin{1});
        end
    end
end

axcodes = {};
for index=1:size(ornt,1)
    axno = ornt(index,1);
    direction = ornt(index,2);
    if isnan(axno)
        axcodes = [axcodes, "None"];
        continue
    end
    axint = int8(round(axno));
    if axint~=axno
        error('ValueError : Non integer axis number %f.', axno);
    elseif direction==1
        axcode = string(labels{axint}{2});
    elseif direction==-1
        axcode = string(labels{axint}{1});
    else
        error('ValueError : Direction should be -1 or 1.');
    end
    axcodes = [axcodes, axcode];
end

end