%%    Overlapper
%     ––––––––––
% Developer: Gregorio Marchesini 
% Date: 9 March 2021

% This script is able to overlap two images taken from the Neptun
% observational campaign. All the images are obtained for HST.

clear all
close all
clc
fprintf('Lading Images....\n')

%% Folders processing and Files Retrival 

files1 = dir('filters');             % add Images Path 
files2 = dir('HST');                 % add Images Path

% withdraw images from the files
files_name = {files1(3:end).name};                                       % The first name is the Strontium fluoride             
addpath(files1(1).folder);                                             % the second is the clear MAMA sensor    
addpath(files2(1).folder);                                             % the second is the clear MAMA sensor    

%    INSERT THE IMAGES HERE 
%–––––––––––––––––––––––––––––––

image_neptune_clear = fitsread('odq408rcq_flt.fits','image',1);  % image darectly from neptune gallery
image_neptune_stf2  = fitsread('odq408raq_flt.fits','image',1);  % image darectly from neptune gallery

% Create two useful copies

image_neptune_lya= image_neptune_clear;
image_neptune_stf2_move  = image_neptune_stf2;


%Grafical object definitiion
HST_image=figure('Position',[0 0 1200 500]);

HST1=subplot(121);hold on
HST2=subplot(122);hold on


HST1.XLim=[350 650];
HST1.YLim=[300 500];
HST1.Title.String={'Clear Filter Image'};
HST1.XLabel.String='pixel';
HST1.YLabel.String='pixel';
colormap(HST1,'jet')
h1=colorbar(HST1);

% this normalises the color bar of the two images
caxis manual
caxis([min(min(image_neptune_clear))    max(max(image_neptune_clear))]);
ylabel(h1,'counts');




HST2.XLim=[350 650];
HST2.YLim=[300 500];
HST2.Title.String={'STF_2 Filter image'};
HST2.XLabel.String='pixel';
HST2.YLabel.String='pixel';
colormap(HST2,'jet')
h2=colorbar(HST2);
ylabel(h2,'counts');

% this normalises the color bar of the two images
caxis manual
caxis([min(min(image_neptune_clear))    max(max(image_neptune_clear))]);


imagesc(HST1,image_neptune_clear)
imagesc(HST2,image_neptune_stf2);




%% Coordinate of Neptune in two different frames 

% I will move the first figure over the second figure

xstf2_centre    = 534;  %pixel
ystf2_centre    = 398;  %pxeli
radius_stf2     = 57 ;  %pixel




xclear_centre     = 491;
yclear_centre     = 380;
radius_clear      = 50 ;

horizontal_move  = xclear_centre-xstf2_centre ;     % movement needed to reovelap the images
vertical_move    = yclear_centre-ystf2_centre ;     % movement needed to reoverlap the images

[row,col]        = size(image_neptune_lya);           % both images are the same size
versor           = sign([horizontal_move,vertical_move]);

horizontal_move  = abs(horizontal_move);              % we are interested in unsig the 
vertical_move    = abs(vertical_move);                % absolute value once we know the sign



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim will be to move the first image (image_neptune_stf2) over the
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
    image_neptune_stf2_move=image_neptune_stf2_move(:,1:col-horizontal_move);  % I cut columns on the right
    image_neptune_stf2_move=image_neptune_stf2_move(1:row-vertical_move,:);    % I cut rows on the lower side

elseif versor == [-1 1]
    image_neptune_stf2_move=image_neptune_stf2_move(:,horizontal_move+1:col);  % I cut the columns on the left
    image_neptune_stf2_move=image_neptune_stf2_move(1:row-vertical_move,:);    % I cut rows on the lower side

elseif versor == [-1 -1]
    image_neptune_stf2_move=image_neptune_stf2_move(:,horizontal_move+1:col);  % I cut the columns on the left
    image_neptune_stf2_move=image_neptune_stf2_move(vertical_move+1:row,:);    % I cut rows on the upper side

elseif versor == [1 -1]
    image_neptune_stf2_move=image_neptune_stf2_move(:,1:col-horizontal_move);  % I cut columns on the right
    image_neptune_stf2_move=image_neptune_stf2_move(vertical_move+1:row,:);    % I cut rows on the upper side
    
end



%% UNIT CONVERSION FROM COUNTS to Raylights

clear_info = fitsinfo('odq408rcq_flt.fits');

% The throughput associated with the two images is different and it is
% obtained in the 'image_processing' script. Hereby the value found in the 
% aforementioned script is used.

scale_factor = 1.1875;                              % clear_image throughput / stf2_image throughput (1600 wavelength)
[t_exp,~]  = KeyFinder(clear_info,'TEXPTIME');      % (s) exposition time 
A_HST      = (2.4/2)^2*pi;                          % [m^2]   HST Mirror aperture
mx         = 0.024;                                 % [arcsec] plate scale value in x direction for a pixel 
my         = 0.024;                                 % [arcsec] plate scale value in y direction for a pixel 
solid      = mx*my*(2*pi/(3600*360))^2              % [steradiants] sterandiant covered by a single pixel
T_lambda   = 0.0387;                                % throughput at 1216 for the clear image
A_eff      = A_HST*T_lambda                         % [m^2]Effective area
counts2ry  = 4*pi*10^-6*(A_eff*10^4)^-2*solid^-1*t_exp^-1; % counts to raylaight converter

% Ry = [photon/s/m^2/cm^2/str]

    %% Overlapping Phase and Subtraction
    
  % in this phase the images are overlapped and subtracted together.
  % the scale factor is used to take into account the fact that the 
  % stf2 image was taken with a filter with different throughput in 
  % comparions with the clear image
    
if     versor == [1 1]
    image_neptune_lya(vertical_move+1:end,horizontal_move+1:end) = ...
    image_neptune_lya(vertical_move+1:end,horizontal_move+1:end) - image_neptune_stf2_move*scale_factor;

elseif versor == [-1 1]   
    image_neptune_lya(vertical_move+1:row,1:col-horizontal_move)  = ...
    image_neptune_lya(vertical_move+1:row,1:col-horizontal_move)  - image_neptune_stf2_move*scale_factor;   

elseif versor == [-1 -1]
    image_neptune_lya(1:row-vertical_move,1:col-horizontal_move)  = ...
    image_neptune_lya(1:row-vertical_move,1:col-horizontal_move)  - image_neptune_stf2_move*scale_factor;

elseif versor == [1 -1]
    image_neptune_lya(1:row-vertical_move,horizontal_move+1:end)  = ...
    image_neptune_lya(1:row-vertical_move,horizontal_move+1:end)  - image_neptune_stf2_move*scale_factor;     
end 
    
%  counts to Raylight Convertion and trimming Phase
%
image_neptune_lya                            = image_neptune_lya*counts2ry;
image_neptune_lya(1:1+vertical_move,:)       = 0;
image_neptune_lya(end-vertical_move:end,:)   = 0;
image_neptune_lya(:,1:horizontal_move)       = 0;
image_neptune_lya(:,end-horizontal_move:end) = 0;

%Grafical object definitiion
comparison=figure('Position',[0 0 1200 500]);
cmp1=subplot(121);hold on
cmp2=subplot(122);hold on


cmp1.XLim          = [350 650];
cmp1.YLim          = [300 500];
cmp1.Title.String  = {'STF_2 and clear after overlapping'};
cmp1.XLabel.String = 'pixel';
cmp1.YLabel.String = 'pixel';

colormap(cmp1,'jet')
h1 =colorbar(cmp1);
ylabel(h1,'Raylight')

% this normalise the color bar of the two images
caxis manual
caxis([min(min(image_neptune_lya))    max(max(image_neptune_lya))]);



cmp2.XLim          = [350 650];
cmp2.YLim          = [300 500];
cmp2.Title.String  = {'STF_2 and clear before corrections '};
cmp2.XLabel.String = 'pixel';
cmp2.YLabel.String = 'pixel';

colormap(cmp2,'jet')
h2= colorbar(cmp2);
ylabel(h2,'Raylight')

% this normalise the color bar of the two images
caxis manual
caxis([min(min(image_neptune_lya))    max(max(image_neptune_lya))]);


imagesc(cmp1,image_neptune_lya)
imagesc(cmp2,(image_neptune_clear-image_neptune_stf2*scale_factor)*counts2ry);

%% Scanning Horizonthal Phase

pixel_range = 1:1024;    % pixel dimension on the image 
range       = 350;        % pixel from the center

x_scanning=sum(image_neptune_lya(yclear_centre-range:yclear_centre+range,xclear_centre-range:xclear_centre+range))*10^-3;   % kR Kiloreylight
y_scanning=sum(image_neptune_lya(yclear_centre-range:yclear_centre+range,xclear_centre-range:xclear_centre+range),2)*10^-3; % kR Kiloreylight


% Graphics

scanner=figure('Position',[0 0 800 400]); 
scanner_ax=axes(scanner);hold on
scatter(cmp1,xclear_centre,yclear_centre,200,'r','filled')
scanner_ax.XGrid='on';
scanner_ax.YGrid='on'
scanner_ax.XLabel.String='Pixel';
scanner_ax.YLabel.String='kR';


plot(xclear_centre-range:xclear_centre+range,x_scanning+0.1,'r','linewidth',3,'DisplayName','raw counts');
plot(yclear_centre-range:yclear_centre+range,y_scanning,'b','linewidth',3,'DisplayName','columns counts');

scatter(xclear_centre,0,200,'filled','HandleVisibility','off')
scatter(yclear_centre,0,200,'filled','HandleVisibility','off')

xline(xclear_centre,'-','x pixel center','LabelVerticalAlignment','bottom','HandleVisibility','off')
xline(yclear_centre,'-','y pixel center','LabelVerticalAlignment','bottom','HandleVisibility','off')
yline(0.1,'-','linewidth',3,'HandleVisibility','off')

legend


% %% Circular Scanning
% 
% for ii = 2:200;
%      mask=circular_kernel(ii,ii-1);
%      % mask=fspecial('disk',ii);
%      brightness(ii)=sum(sum(mask.*image_neptune_lya(yclear_centre-ii:yclear_centre+ii,xclear_centre-ii:xclear_centre+ii)*10^-3));
% end
% 
% figure()
% plot(brightness)






%%
fprintf('Finished of the process....\nEnjoy \n')








