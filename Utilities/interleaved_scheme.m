%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that returns the indexes of the slices to be acquired in an   %
%  interleaved manner depending on the number of slices.                  %
%                                                                         %
%         interleavedSlices_index = interleaved_scheme(NbSlices);         %
%                                                                         %
%  input:   - NbSlices: number of slices                                  %
%                                                                         %
%  output:  - interleavedSlices_index: list of indexes of the slices      %
%                                      corresponding to an interleaved    %
%                                      acquisition scheme                 %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2021-04-21                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function interleavedSlices_index = interleaved_scheme(NbSlices)

% Input check
if nargin < 1
    error('Missing input(s).');
elseif nargin > 1
    error('Too many inputs.');
end

% Implement interleaved slice acquisition scheme
if mod(NbSlices,2)==0
	interleavedSlices_index = [2:2:NbSlices, 1:2:NbSlices-1];
else
    interleavedSlices_index = [2:2:NbSlices-1, 1:2:NbSlices];
end

% Display message for debugging
sprintf('Slices will be acquired in the following order:')
fprintf(' %d\n', interleavedSlices_index(:))

end