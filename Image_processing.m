% This script is intended to retrive flux information from the photon
% counter device MAMA which is used by STIS in order to obtain images of
% the planet Neptune. The script retrives the images from a specific 
% fits file where of the '_flt' extension (flat filed science). The images 
%                                                            
% given observations are in counts (number of photons hitting the sensor)
% we would like to be able to convert this count in flux units (Jensky) or
% in Reyleight.

% The problem consist in stating valuable assumptions in order to connect
% the counts to the wavelenght of the photons hitting the camera. In
% particular the energy of the photon hitting the sensor is directly link
% to his energy and the total power received from the MAMA sensor can be
% estimated only if the wavelength of the incident light is known. A seires
% of filters can be applied in order to select the desired wavelength to 
% obsereve. However not all the filters have the same efficiency. With the
% term efficinecy we refer to the fact that not all ther light enetering
% the telescopoe ca effectively reach the sensor, due to losses inside the
% optical system and the filter itself. This can significantly affect the 
% estimated flux measured through a count to energy conversion.
% In the specific case of study, the Lyman-alpha emission (1216 A)
% is the target wavelength.

% to be continued

% Strontioum-Fluoride filter (StF2) throughput info:
% https://hst-docs.stsci.edu/stisihb/chapter-14-imaging-reference-material/14-5-fuv-mama/f25srf2-fuv-mama-longpass
clc
clear all
close all


%% Initial Parameters
cal_wavelength=1600;    % [Å] Amstrong  Throughput calibration wavelength
Lya_wavevlenght=1216;   % [Å] Amstrong  Target wavelngth is the lyman alpha


%% Filter Parameters
files1=dir('filters');
files2=dir('HST');                 % the first two names are always '..' and '.' 
                                   % because they makes you come back to the directory
                                   % remeber to eliminate the if you want
                                   % to use them
                                                           
disp('-------------------------------')
disp('LIST OF FILES')
disp('-------------------------------')
fprintf('%s\n',files2(3:end).name)
disp('-------------------------------')

files_name={files1(3:end).name};   % The first name is the Strontium fluoride             
addpath(files1(1).folder);         % the second is the clear MAMA sensor    
addpath(files2(1).folder);         % the second is the clear MAMA sensor    

filter1=dlmread(files_name{1});
filter2=dlmread(files_name{2});

%% Plots for the Figure

image_fig=figure();

% Clear MAMA Detector

subplot(4,2,1)
hold on
x1=filter1(:,1);                             % wavelength
y1=filter1(:,2);                             % percetage throughput
plot(x1,y1,'r','linewidth',4);
target_thp_Cl=y1(x1==cal_wavelength);        % target through put at target_wavelength (Å)

scatter(cal_wavelength,y1(x1==cal_wavelength),100,'b','linewidth',8)
string_format=sprintf('Throughput: \n %.3f',y1(x1==cal_wavelength));
text(cal_wavelength+100,y1(x1==cal_wavelength)+0.005,string_format)
grid on
xlabel('Wavelength (Å)')
ylabel('Throughput (%)')
title(sprintf(files_name{1}),'Interpreter', 'none')

% Strontium Fluoride Filter

subplot(4,2,2)
hold on
x2=filter2(:,1);                               % wavelength (amstrong)
y2=filter2(:,2);                               % percetage throughput %100)
plot(x2,y2,'b','linewidth',4);
target_thp_StF2=y2(x2==cal_wavelength);        % target through put at target_wavelength (Å)

scatter(cal_wavelength,y2(x2==cal_wavelength),100,'b','linewidth',8)
string_format=sprintf('Throughput : \n %.3f',y2(x2==cal_wavelength));
text(cal_wavelength+100,y2(x2==cal_wavelength)+0.005,string_format)
grid on
xlabel('Wavelength (Å)')
ylabel('Throughput (%)')
title(sprintf(files_name{2}),'Interpreter', 'none')


%% Parameters Identification

%Throughput Ratio

scaling_factor=target_thp_Cl/target_thp_StF2;  % this is the normalization factor for
                                               % our image. This is because the
                                               % througputs are different for the two
                                               % images that must be
                                               % subtracted

                                               
             
disp('Loading Requested Images ...')                                   
%% FITS file reading

fits_info_stf2=fitsinfo('odq405acq_flt.fits');                     % info file
filter_name_stf2=KeyFinder(fits_info_stf2,'FILTER');               % read the filter used
filter_exptime_stf2=KeyFinder(fits_info_stf2,'TEXPTIME');          % read the exposure time [s]
image_neptune_stf2=fitsread('odq405acq_flt.fits','image',1);       % image reading

fits_info_clear=fitsinfo('odq405aeq_flt.fits');                    % info file
filter_name_clear=KeyFinder(fits_info_clear,'FILTER');             % read the filter used
filter_exptime_clear=KeyFinder(fits_info_clear,'TEXPTIME');        % read the exposure time [s]
image_neptune_clear=fitsread('odq405aeq_flt.fits','image',1);      % image reading

%% Throughput Normalization


image_neptune_clear=image_neptune_clear*scaling_factor;  % This operation equates
                                                         % the throughputs 
                                                         % of the two images so
                                                         % that the final result
                                                         % takes into account for
                                                         % the lower throughput of
                                                         % the filtered image 

LYalpha_image=image_neptune_clear-image_neptune_stf2;
% This image has (ideally) only lyman alpha radiation after the subtraction
% The first image contains all the wavelengths. The second contains all the
% wavelengths except for the Lyman Alpha.

%% Plot The resuls from the subtraction

subplot(4,2,[3 4 5 6 7 8])
imagesc(LYalpha_image)
set(image_fig,'Position',[0,0,600,800]);
title(sprintf('%s - %s',filter_name_clear,filter_name_stf2),'Interpreter','none')

figure()

mesh(LYalpha_image)
xlabel('pixel')
ylabel('pixel')
title('Lyman_\alpha image')

figure()

mesh(image_neptune_clear)
xlabel('pixel')
ylabel('pixel')
title('Clear Image')

figure()

mesh(image_neptune_stf2)
xlabel('pixel')
ylabel('pixel')
title('StF_2 Filter')


% Some comments on the 3D results:
% from the 3D image of the stf2 filter it seems that i has way less counts
% distributed on the image. Is this a sign that most of the lyght is in
% lyman alpha? because to me the two filters should be equally populated

disp('Loading Finished ...') 


%% Conversion Between Photon Count and Energy

h_p       = 6.62607004e-34             % [m^3 kg s^-1] plank constant
c         = 3e8;                       % [m s^-1]      speed of light
A2m       = 10^-10;                    % Å to m conversion
E_lyman   = h_p*c/Lya_wavevlenght;     % Energy of one photon in Lyman Alpha wavelengh [J]
pow_lyman = E_lyman/1

% %% Converted image
% 
% figure()
% Lyman_energy=imagesc(LYalpha_image)*E_lyman;


  