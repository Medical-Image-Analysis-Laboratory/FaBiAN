%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that saves the T2w simulation of the fetal brain along with   %
%  its label map according to the gestational age of the fetus and to     %
%  the acquisition plane                                                  %
%                                                                         %
%                 run_id = name_run(orientation, shift_mm)                %
%                                                                         %
%  inputs:  - orientation: strict acquisition plane (axial, coronal or    %
%                          sagittal)                                      %
%           - shift_mm: displacement (in mm) of the slice slab in the     %
%                       slice thickness direction                         %
%                                                                         %
%                                                                         %
%  Hélène Lajous, 2022-12-13                                              %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_id = name_run(orientation, shift_mm)

% Input check
if nargin < 2
    error('Missing input(s).');
elseif nargin > 2
    error('Too many inputs.');
end

if orientation==1
    if shift_mm==0
        run_id = 1;
    elseif shift_mm<0
        run_id = 2;
    elseif shift_mm>0
        run_id = 3;
    end
elseif orientation==2
    if shift_mm==0
        run_id = 4;
    elseif shift_mm<0
        run_id = 5;
    elseif shift_mm>0
        run_id = 6;
    end
elseif orientation==3
    if shift_mm==0
        run_id = 7;
    elseif shift_mm<0
        run_id = 8;
    elseif shift_mm>0
        run_id = 9;
    end
end

end