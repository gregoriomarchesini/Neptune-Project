function string=s2h(ts)

% Developer   : Gregorio Marchesini
% Date        : 02/2021
% Description : Converts an interval of time from seconds to
%               hours,minutes,second format

% Description
% –––––––––––

% This function converts seconds into days-hour-minutes-seconds format.
% The output value is a string

% INPUTS
% ––––––

%   Name     Data Type     Description 

%   ts       scalar        time in seconds

% OUTPUT
% ––––––––

%   Name     Data Type     Description 

%   string     string      Returned fromatted time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FUNCTION
  ts     = floor(ts);
  h      = floor(ts/60/60);
  m      = floor((ts-h*60*60)/60);
  s      = ts-((h*60*60)+(m*60));
  string = [num2str(h),'h-',num2str(m),'m-',num2str(s),'s'];
end