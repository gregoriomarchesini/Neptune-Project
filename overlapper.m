%%    Overlapper
%     ––––––––––
% Developer: Gregorio Marchesini 
% Date: 9 March 2021

% This script is able to overlap two images taken from the Neptun
% observational campaign. All the images are obtained for HST.

clear all
close all
clc
fprintf('Loading Images....\n')

%% Folders processing and Files Retrival 

files1 = dir('filters');             % add Images Path 
files2 = dir('HST');                 % add Images Path

% withdraw images from the files

files_name = {files1(3:end).name};                                       % The first name is the Strontium fluoride             
addpath(files1(1).folder);                                             % the second is the clear MAMA sensor    
addpath(files2(1).folder);                                             % the second is the clear MAMA sensor    

%    INSERT THE IMAGES HERE 
%–––––––––––––––––––––––––––––––

filename_clear      = 'odq406twq_flt.fits';
filename_srf2       = 'odq406tuq_flt.fits';


image_neptune_clear = fitsread(filename_clear,'image',1);  % image darectly from neptune gallery
image_neptune_srf2  = fitsread(filename_srf2,'image',1);   % image darectly from neptune gallery

clear_info          = fitsinfo(filename_clear);            % input informations
srf2_info           = fitsinfo(filename_srf2);             % input informations

% exposure time information

[t_exp_clear,~]    = KeyFinder(clear_info,'TEXPTIME');      % (s) exposure time clear image
[t_exp_srf2,~]     = KeyFinder(srf2_info,'TEXPTIME') ;      % (s) exposure time srf2 image   

% NORMILIZE FOR THE EXPOSITION TIME : from counts to counts/s

image_neptune_clear = image_neptune_clear/t_exp_clear;   % [counts/s]
image_neptune_srf2  = image_neptune_srf2/t_exp_srf2;     % [counts/s]

image_neptune_clear(image_neptune_clear < 0) = 0;        % Eliminate the negative counts
image_neptune_srf2(image_neptune_srf2 < 0)   = 0;        % Eliminate the negative counts
 
% Create two useful copies

image_neptune_lya        = image_neptune_clear; 
image_neptune_srf2_move  = image_neptune_srf2; 


%Grafical object definitiion
HST_image=figure('Position',[0 0 1200 800]);

% genearl title
sgtitle(filename_clear,'interpreter','none')

% subplot axes

HST1=subplot(221);hold on
HST2=subplot(222);hold on


HST1.XLim=[350 650];
HST1.YLim=[300 500];
HST1.Title.String={'Clear Filter Image'};
HST1.XLabel.String='pixel';
HST1.YLabel.String='pixel';
h1=colorbar(HST1);
colormap(HST1,'hot')
% this normalises the color bar of the two images
caxis manual
caxis([min(min(image_neptune_clear))    max(max(image_neptune_clear))]);
ylabel(h1,'counts/s');


HST2.XLim=[350 650];
HST2.YLim=[300 500];
HST2.Title.String={'STF_2 Filter image'};
HST2.XLabel.String='pixel';
HST2.YLabel.String='pixel';
h2=colorbar(HST2);
colormap(HST2,'hot')
ylabel(h2,'counts/s');

% this normalises the color bar of the two images
caxis manual
caxis([min(min(image_neptune_clear))    max(max(image_neptune_clear))]);

imagesc(HST1,image_neptune_clear);  % counts/s
imagesc(HST2,image_neptune_srf2);   % counts/s


%% Coordinate of Neptune in two different frames 

% I will move the first figure over the second figure

xsrf2_centre      = 535 ;  %pixel
ysrf2_centre      = 399;  %pxeli
radius_srf2       = 55 ;  %pixel


xclear_centre     = 493 ;
yclear_centre     = 379;
radius_clear      = 61 ;

horizontal_move  = xclear_centre-xsrf2_centre ;       % movement needed to reovelap the images
vertical_move    = yclear_centre-ysrf2_centre ;       % movement needed to reoverlap the images

[row,col]        = size(image_neptune_lya);           % both images are the same size
versor           = sign([horizontal_move,vertical_move]);

horizontal_move  = abs(horizontal_move);              % we are interested in unsig the 
vertical_move    = abs(vertical_move);                % absolute value once we know the sign



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim will be to move the first image (image_neptune_srf2) over the
% second one in order to achieve perfect overlapping of the planet neptune
% in the two images. This wioll come with a lost of information at the
% edges of the image, but not a lot of bits should be lost in the process.
%The edges of the moving images will be just trimmed a bit but the planet
%remains our interst
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trimming phase

% two possible movements are possible for eax axis. I can move the image 
% up or down and left or right. The sign of horizonatal_move and
% vertical_move will tell me what is going on. 

% horizonatal_move > 0  ---->  I will need to move to the right
% horizonatal_move < 0  ---->  I will need to move to the left
% vertical_move    > 0  ---->  I will need to move to the up
% vertical_move    < 0  ---->  I will need to move to the down


% the axes of refernce are the clear image matrix axes of reference

%      1 2 3 4 5 6 7 .....
%      1
%      2
%      3
%      4
%      5
%      :
%      :

% for each case I will have to trimm in a different way. The moving image
% is always the one who is trimmed.
% NOTE: the case in which the movemnet is zero leaves the axis unaffected



if     versor == [1 1]
    image_neptune_srf2_move=image_neptune_srf2_move(:,1:col-horizontal_move);  % I cut columns on the right
    image_neptune_srf2_move=image_neptune_srf2_move(1:row-vertical_move,:);    % I cut rows on the lower side

elseif versor == [-1 1]
    image_neptune_srf2_move=image_neptune_srf2_move(:,horizontal_move+1:col);  % I cut the columns on the left
    image_neptune_srf2_move=image_neptune_srf2_move(1:row-vertical_move,:);    % I cut rows on the lower side

elseif versor == [-1 -1]
    image_neptune_srf2_move=image_neptune_srf2_move(:,horizontal_move+1:col);  % I cut the columns on the left
    image_neptune_srf2_move=image_neptune_srf2_move(vertical_move+1:row,:);    % I cut rows on the upper side

elseif versor == [1 -1]
    image_neptune_srf2_move=image_neptune_srf2_move(:,1:col-horizontal_move);  % I cut columns on the right
    image_neptune_srf2_move=image_neptune_srf2_move(vertical_move+1:row,:);    % I cut rows on the upper side
    
end



%% UNIT CONVERSION FROM COUNTS to Raylights

% The throughput associated with the two images is different and it is
% obtained in the 'image_processing' script. Hereby the value found in the 
% aforementioned script is used.

scale_factor = 1.1875;                                    % clear_image throughput / srf2_image throughput (1600 wavelength) % Remmeber we assume only two wavelengths are present in the clear image (2600 and 1216) while only one is present in the srf2 (1600)
A_HST        = (2.4/2)^2*pi;                              % [m^2]   HST Mirror aperture
mx           = 0.0246;                                    % [arcsec] plate scale value in x direction for a pixel 
my           = 0.0246;                                    % [arcsec] plate scale value in y direction for a pixel 
solid        = mx*my*(2*pi/(3600*360))^2;                 % [steradiants] sterandiant covered by a single pixel
T_lambda     = 0.0387;                                    % throughput at 1216 for the clear image
A_eff        = A_HST*T_lambda;                            % [m^2]Effective area
counts2ry    = 4*pi*10^-6*(A_eff*10^4)^-1*solid^-1*10^-3; % counts/s to kilo raylaight converter [Brightness]

% REMEMBER: Ry = [photon/s/m^2/cm^2/str] [BRIGHTNESS]

  %% Overlapping Phase and Subtraction
    
  % in this phase the images are overlapped and subtracted together.
  % the scale factor is used to take into account the fact that the 
  % srf2 image was taken with a filter with different throughput in 
  % comparions with the clear image
    
if     versor == [1 1]
    image_neptune_lya(vertical_move+1:end,horizontal_move+1:end) = ...
    image_neptune_lya(vertical_move+1:end,horizontal_move+1:end) - image_neptune_srf2_move*scale_factor;

elseif versor == [-1 1]   
    image_neptune_lya(vertical_move+1:row,1:col-horizontal_move)  = ...
    image_neptune_lya(vertical_move+1:row,1:col-horizontal_move)  - image_neptune_srf2_move*scale_factor;   

elseif versor == [-1 -1]
    image_neptune_lya(1:row-vertical_move,1:col-horizontal_move)  = ...
    image_neptune_lya(1:row-vertical_move,1:col-horizontal_move)  - image_neptune_srf2_move*scale_factor;

elseif versor == [1 -1]
    image_neptune_lya(1:row-vertical_move,horizontal_move+1:end)  = ...
    image_neptune_lya(1:row-vertical_move,horizontal_move+1:end)  - image_neptune_srf2_move*scale_factor;     
end 
    
%  counts to kiloRaylight Convertion and trimming Phase

image_neptune_lya                            = image_neptune_lya*counts2ry;
image_neptune_lya(1:1+vertical_move,:)       = 0;
image_neptune_lya(end-vertical_move:end,:)   = 0;
image_neptune_lya(:,1:horizontal_move)       = 0;
image_neptune_lya(:,end-horizontal_move:end) = 0;

%Grafical object definitiion


cmp1=subplot(223);hold on
cmp2=subplot(224);hold on


cmp1.XLim          = [350 650];
cmp1.YLim          = [300 500];
cmp1.Title.String  = {'STF_2 and clear after overlapping and subtraction','Only Lyman-\alpha remaining'};
cmp1.XLabel.String = 'pixel';
cmp1.YLabel.String = 'pixel';

colormap(cmp1,'hot')
h1 =colorbar(cmp1);
ylabel(h1,'Brightness [kR]')

% this normalises the color bar of the two images
caxis manual
caxis([min(min(image_neptune_lya))    max(max(image_neptune_lya))]);



cmp2.XLim          = [350 650];
cmp2.YLim          = [300 500];
cmp2.Title.String  = {'STF_2 and clear before corrections and subtraction '};
cmp2.XLabel.String = 'pixel';
cmp2.YLabel.String = 'pixel';

colormap(cmp2,'hot')
h2= colorbar(cmp2);
ylabel(h2,'Brightness [kR]')

% this normalise the color bar of the two images

caxis manual
caxis([min(min(image_neptune_lya))    max(max(image_neptune_lya))]);

imagesc(cmp1,image_neptune_lya)
imagesc(cmp2,(image_neptune_clear-image_neptune_srf2*scale_factor)*counts2ry);



%% Scanning Horizonthal Phase

pixel_range     = 1:1024;     % pixel dimension on the image 
xrange          = 200;        % pixel from the center
yrange          = 50;
x_scanning_mean = mean(image_neptune_lya(abs(yclear_centre-yrange):abs(yclear_centre+yrange),abs(xclear_centre-xrange):abs(xclear_centre+xrange)));   % kR Kiloreylight summation column-wise
y_scanning_mean = mean(image_neptune_lya(abs(yclear_centre-yrange):abs(yclear_centre+yrange),abs(xclear_centre-xrange):abs(xclear_centre+xrange)),2); % kR Kiloreylight summation raw-wise


% Graphics
rectangle(cmp1,'Position',[xclear_centre-200 yclear_centre-50 400 100])
scanner=figure('Position',[0 0 800 400]); 
scanner_ax=axes(scanner);hold on
scatter(cmp1,xclear_centre,yclear_centre,100,'k','filled')
scanner_ax.XGrid='on';
scanner_ax.YGrid='on';
scanner_ax.XLabel.String='Pixel';
scanner_ax.YLabel.String='mean Brightness [kR]';

% NOTE THE FIRST PLOT IS RISED OF 1kR for good separation of the plots

offset = 1;       %  [KR]
plot(scanner_ax,xclear_centre-xrange:xclear_centre+xrange,x_scanning_mean+offset,'r','linewidth',3,'DisplayName','mean along the columns');
plot(scanner_ax,yclear_centre-yrange:yclear_centre+yrange,y_scanning_mean,'b','linewidth',3,'DisplayName','mean along the raws');
scatter(scanner_ax,xclear_centre,0,200,'filled','HandleVisibility','off')
scatter(scanner_ax,yclear_centre,0,200,'filled','HandleVisibility','off')

xline(scanner_ax,xclear_centre,'-','x pixel center','LabelVerticalAlignment','bottom','HandleVisibility','off')
xline(scanner_ax,yclear_centre,'-','y pixel center','LabelVerticalAlignment','bottom','HandleVisibility','off')
yline(scanner_ax,offset,'-','offset of 1 kR on the red plot','linewidth',3,'HandleVisibility','off')

legend

%%
% Scanning for average Lyman_alpha radiation

R_neptune   = radius_clear;                      % how much further is the scanning going from the center of the planet
inner_limit = 1;                                 % inner limit
counter     = 1;
inter       = 2;
range01= inner_limit:inter:R_neptune;
range12= R_neptune+1:inter:2*R_neptune;
range23= 2*R_neptune+1:inter:5*R_neptune;

for ii = range01
  Area_disk                        = (ii^2)*pi;    
  planet_disk_kernel               = circular_kernel(ii,0);                  % disk kernel
  brightness_planet(counter)       = sum(sum(planet_disk_kernel.*image_neptune_lya(yclear_centre-ii:yclear_centre+ii,xclear_centre-ii:xclear_centre+ii)))/Area_disk;
  counter                          = counter+1;
end
  counter     = 1;
for ii = range12
   Area_anulus_12                  = ((ii)^2-(range12(1)-1)^2)*pi;     % anulus area covered by the kernel from R_nep to 2xR_nep
   bkg_anulus_kernel               = circular_kernel(ii,range12(1)-1);
   brightness_bkg12(counter)       = sum(sum(bkg_anulus_kernel.*image_neptune_lya(yclear_centre-(ii):yclear_centre+(ii),xclear_centre-(ii):xclear_centre+(ii))))/Area_anulus_12;     % background mean value over an Anulus around the planet
   counter                         = counter+1;
end
counter     = 1;

for ii = range23
   Area_anulus_23                 = ((ii)^2-(range23(1)-1)^2)*pi;     % anulus area covered by the kernel from R_nep to 2xR_nep
   bkg_anulus_kernel              = circular_kernel(ii,range23(1)-1);
   brightness_bkg23(counter)      = sum(sum(bkg_anulus_kernel.*image_neptune_lya(yclear_centre-(ii):yclear_centre+(ii),xclear_centre-(ii):xclear_centre+(ii))))/Area_anulus_23;     % background mean value over an Anulus around the planet
   counter                        = counter+1;
end




%brightness density over the disk

figure('Position',[50 50 1200 400])
miny = min(brightness_planet);
maxy = max(brightness_planet);
ax1=subplot(131);hold on
ax1.YLabel.String='Average Brightness over the Disk from 0 to 1 R_{N} [kR]';
ax1.XLabel.String='Pixel offset from the center of the planet ';
ax1.YLim = [0 maxy];
grid on

plot(ax1,range01,brightness_planet,'k','linewidth',2)
xline(ax1,R_neptune,'r--','R_{nep}')

ax2=subplot(132);hold on
ax2.YLabel.String='Average Brightness over the Anulus from 1 R_{N} to 2 R_{N} [kR]';
ax2.XLabel.String='Pixel offset from the center of the planet ';
ax2.YLim = [0 maxy];

plot(ax2,range12,brightness_bkg12,'k','linewidth',2)
xline(ax2,R_neptune,'r--','R_{nep}')
xline(ax2,2*R_neptune,'r--','2xR_{nep}')
grid on

ax3=subplot(133);hold on
ax3.YLabel.String='Average Brightness over the Anulus from 2 R_{N} to 3 R_{N} [kR]';
ax3.XLabel.String='Pixel offset from the center of the planet ';
ax3.YLim = [0 maxy];

plot(ax3,range23,brightness_bkg23,'k','linewidth',2)
xline(ax3,2*R_neptune,'r--','2xR_{nep}')
xline(ax3,ceil(range23(end)/R_neptune)*R_neptune,'r--',sprintf('%ixR_{nep}',ceil(range23(end)/R_neptune)))
grid on




%% Polynomial fitting
figure
ax=axes();
holed_image=image_neptune_lya;
holed_image(yclear_centre-radius_clear:yclear_centre+radius_clear,xclear_centre-radius_clear:xclear_centre+radius_clear)=0;
imagesc(holed_image)
colormap(ax,'jet')
h2 = colorbar(ax);
ylabel(h2,'Brightness [kR]')

fprintf('Finish of the process....\nEnjoy \n')


%% Plinomial Curve fitting 





