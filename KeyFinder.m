function [result,description]=KeyFinder(fits_info_obj,desired_key)

% Developer: Gregorio Marchesini 
% Date: 1 March 2021


% This function was create in order to easily obtain Keyword values from
% the Primary Header of a fits file. Using arrayfun it is possible to  give multiple
% words input to this function if necessary. Otherwise you can simply
% insert this ina loop function and obtain all the desired words

%  INPUTS
%  ––––––
% 
% 
%  Variable           Description                   Data Type
% ––––––––––          ––––––––––––                  ––––––––––
%
% fits_info_obj     object obtained reading         fits-object
%                   a fits file with fitsinfo
%
% desired_key       Keyword to reseacrh in the      string
%                   primary header
%
%-------------------------------------------------------------------------
%
%  OUTPUTS
%  –––––––
%
%  Variable           Description                   Data Type
% ––––––––––          ––––––––––––                  ––––––––––
%
%  result             String related to Keyword     string
%                     selected
% 
%
%
%
% NOTE: A complete list of the possible headers can be found in the HST
% Doucumentation. The word must be inserted correctly


% extract cell array of all the the Keywords
table_key    = fits_info_obj.PrimaryData.Keywords; 
[r,~]        = size(table_key);
result       = NaN;
descpription = NaN;
flag         = 0;


for i=1:r
    if strcmp(desired_key,table_key{i,1})
        result      = table_key{i,2};
        description = table_key{i,3};
        flag=1;
    end  
end

if flag==0
    fprintf('Sorry, no results for yout search\n')
end

end
