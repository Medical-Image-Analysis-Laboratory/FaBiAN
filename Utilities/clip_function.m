%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function that transforms the advantage/disadvantage values to          %
%  new controlled ones through a sigmoid                                  %
%                                                                         % 
%                  NewMap = clip_function(      Map, ...                  %
%                                         ClipValue, ...                  %
%                                                GA)                      %
%                                                                         %
%  inputs:  - Map: advantage/disadvantage values to tune T1/T2 maps       %
%           - ClipValue: control value set not to deviate T1 and T2       %
%                        values more than a specific percentage of the    %
%                        corresponding reference T1, respectively T2      %
%                        value. The 'adapt' option takes the clipping     %
%                        value from a list derived from analyzing the     %
%                        contrast in T2-weighted images of a subject      %
%                        from a normative spatiotemporal MRI atlas of     %
%                        the fetal brain of equivalent gestational age    %
%                        as the subject being simulated.                  %
%           - GA: gestational age of the fetus (in weeks)                 %
%                                                                         %
%                                                                         %
%  outputs: - NewMap: map of controlled advantage/disadvantage values     %
%                                                                         %
%                                                                         %
%  le Boeuf Andrés, 2022-04-22                                            %
%  andres.le.boeuf@estudiantat.upc.edu                                    %
%  Modified by Hélène Lajous, 2023-02-15                                  %
%  helene.lajous@unil.ch                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NewMap = clip_function(Map, ClipValue, GA)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs!');
end

%Clipping value selection
if isequal(ClipValue, 'adapt')  %adapt clipping value to the GA of the subject
    % The following clipping values were determined iteratively based on
    % the histogram distribution of subjects from a normative
    % spatiotemporal MRI atlas (STA) of the fetal brain (Gholipour et al.,
    % Scientific Reports, 2017)
    ClipValues_STA = [0.2, 0.23, 0.26, 0.32, 0.36, 0.38, 0.37, 0.39, 0.39, 0.39, 0.36, 0.37, 0.43, 0.4, 0.37, 0.33, 0.32, 0.26];
    GA_STA = 21:38; %range of gestational age (GA, in weeks) in STA
    % Interpolate clipping values for GA given in a decimal form between
    % 20.0 and 38.0 weeks by shape-preserving piecewise cubic interpolation
    ClipValuesInterpolated = interp1(GA_STA, ClipValues_STA, 20.0:0.1:38.0, 'pchip');
    ClipValueReward = ClipValuesInterpolated(round((GA - 20) * 10 + 1));
elseif (ClipValue>=0 && ClipValue<=1) %if want to specify explicitely the clipping value to use
    ClipValueReward = ClipValue;
else
    error("Wrong clip value.")
end
fprintf("Current clipping value: %f\n", ClipValueReward)

% Apply the clipping value to the advantage/disadvantage values
Map = ClipValueReward * (2 ./ (1 + exp((-2 / ClipValueReward) * Map)) - 1);   %modulate T1/T2 values by sigmoid of parameter clipping value
NewMap = Map + 1;  %because want to implement a variation of mean ref T1 and T2 values (Ã©ventuellement Ã  sortir de cette fonction et passer dans WM_maturation.m)
    
end