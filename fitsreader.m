

function [FITSINFO,FITSDATASET]=fitsreader(FITSDIR,name,varargin)

%% FUNCTION fitsReading
%  version/date : version 01, 131203
%  author(s)    : Fabrizio-M. Musacchio, IGM Cologne
%                 (musacchio@geo.uni-koeln.de.de)
%  modified by Lorenz Roth (lorenzr@kth.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update : This version was last edited by Gregorio Marchesini (gremar@kth.se)
% Date   : 02/2021

%% DESCRIPTION
%
%  INPUT:   FITSDIR         directory of the observation
%           name     Name fits-file.
%  
%  NOTES : the path must be correctly inserted in the function. The
%          observation can be called by name. No extension is needed
%          like '_flt' since this is the only one cupported so far 
%          Pag 31 DATA HANDBOOK STIS (Version April 2019)
%          for further specifications on the possible suffix
%
%
%
%  OUTPUT:  FITSDATASET     1024x1024 FITS-Image + additional Information
%  (header) about dispersion, central wavelength, etc.
%
%  Example:
           %DATADIR     = './data_HST/';
           %OBSERVATION = 'ocai02040';           % root name of file

%           [F] = fitsReading(DATADIR, HSTcampaign, OBSERVATION, VISOBStmp);
%
%           --> F.RAWCOUNTS is then 1024x1024px raw image
%               F.EXPTIME is Exposure Time of image
%               
          
% NOTE: fitsinfo gives you only the info in the extensions 
%       fits read is able to read the content of the object

%% Input Parser 

% Defoult Parameters
printer={'no','yes'};

% Parsing phase
p=inputParser();
addParameter(p,'plots_in','no',@(x)mustBeMember(x,printer)); %parameter plots_in to decide if you want to plot or not
a=parse(p,varargin{:});


% Read FITS-file:

FITSINFO                 = fitsinfo(FITSDIR);
FITSDATASET.RAWCOUNTS    = fitsread(FITSDIR,'image',1);

[name_target,~] = KeyFinder(FITSINFO,'TARGNAME');
[exptime,~]     = KeyFinder(FITSINFO,'TEXPTIME');
[date_obs,~]    = KeyFinder(FITSINFO,'TDATEOBS');
[Filter,~]      = KeyFinder(FITSINFO,'FILTER');
[FOV,~]         = KeyFinder(FITSINFO,'APER_FOV');
[pixel_arc,~]   = KeyFinder(FITSINFO,'PLATESC');



% Read specific parameters from header information
fprintf('Dataset Description : %s \n',name)
disp('--------------------------------------------------------')
fprintf('Target Name                      : %s \n',name_target)
fprintf('Total exposure time              : %s \n',s2h(exptime));         % Exposure time [s]
fprintf('Starting date of the observation : %s\n',date_obs);  % time start exposure
fprintf('Filter Used                      : %s\n',Filter);                   % Filter used in the observation
fprintf('Filed Of View                    : %s\n',FOV)
fprintf('pixel resolution                 : %s pixel/arcsec\n',pixel_arc)
disp('--------------------------------------------------------')
fprintf('\n\n\n\n')

%% Flux Units Conversion  (pag 174 HST DATA HandBOOK version 2019)
photflam=FITSINFO.PrimaryData.Keywords{100,2};    %[erg cm^-2 sec^-1 Å^-1] sensitivity of the sensor
exposure_time=FITSINFO.PrimaryData.Keywords{41,2}; %[s] total exposure time
flux_conv=photflam/exposure_time;


%% PLOTS the image 

%% Figure set up
if strcmp(p.Results.plots_in,'yes')
    figure()
   
    set(gca,'ydir','nor')  % y-Achse soll von unten nach oben verlaufen
    s=imagesc( 1:1024,...
        1:1024,...
        FITSDATASET.RAWCOUNTS*flux_conv);
    
    h=colorbar;
    title(h,'erg cm-^2 s-^1 Å-^1','FontSize',16);
    axis equal
    
    xlim([1 1024]);
    ylim([1 1024]);
    xlabel('pixel');
    ylabel('pixel');
end
end
