%% Overlapper
% This script can overlap the HST images in order to obtain the correct
% datas in terms of overlapping imges after a proper cutting is done
% the coordinate for the center of the two images is obtained bt using 
% the script image_gridder

close all
clc
clf
clear all


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

% I create two useful copies

image_neptune_overlapped = image_neptune_clear;
image_neptune_stf2_move  = image_neptune_stf2;





%% Coordinate of Neptune in two different frames 

% I will move the first figure over the second figure

xclear_centre  = 491;
yclear_centre  = 380

xstf2_centre  = 534;
ystf2_centre  = 398;


horizontal_move  =   xclear_centre-xstf2_centre ;     % movement needed to reovelap the images
vertical_move    =   yclear_centre-ystf2_centre ;     % movement needed to reoverlap the images

[row,col] = size(image_neptune_overlapped);           % both images are the same size
versor    = sign([horizontal_move,vertical_move]);

horizontal_move = abs(horizontal_move);              % we are interested in unsig the 
vertical_move   = abs(vertical_move);                % absolute value once we know the sign





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
% NOTE: the case in which the movemnet is zero leav the axis unaffected



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

    %% Overlapping Phase
    
    % FOR NOW WE ARE TAKING THE SUM OF THE OVERLAPPED IMAGES BUT BE AWARE
    % THAT THE DIFFERNCE SHOULD BE TAKEN IN THE ANALYSIS INSTEAD
    
if     versor == [1 1]
    image_neptune_overlapped(vertical_move+1:end,horizontal_move+1:end) = ...
    image_neptune_overlapped(vertical_move+1:end,horizontal_move+1:end) + image_neptune_stf2_move;

elseif versor == [-1 1]   
    image_neptune_overlapped(vertical_move+1:row,1:col-horizontal_move)  = ...
    image_neptune_overlapped(vertical_move+1:row,1:col-horizontal_move)  + image_neptune_stf2_move;   

elseif versor == [-1 -1]
    image_neptune_overlapped(1:row-vertical_move,1:col-horizontal_move)  = ...
    image_neptune_overlapped(1:row-vertical_move,1:col-horizontal_move)  + image_neptune_stf2_move;

elseif versor == [1 -1]
    image_neptune_overlapped(1:row-vertical_move,horizontal_move+1:end)  = ...
    image_neptune_overlapped(1:row-vertical_move,horizontal_move+1:end)  + image_neptune_stf2_move;     
end 
    
%Grafical object definitiion
comparison=figure('Position',[0 0 1200 500])

cmp1=subplot(131);hold on
cmp2=subplot(132);hold on
cmp3=subplot(133);hold on

cmp1.XLim=[350 650];
cmp1.YLim=[300 500];
cmp1.Title.String={'subtraction overlapped'};
cmp2.XLim=[350 650];
cmp2.YLim=[300 500];
cmp2.Title.String={'no overlapping'};
cmp2.XLim=[350 650];
cmp2.YLim=[300 500];
cmp2.Title.String={'summation overlapping'};

imagesc(cmp1,image_neptune_overlapped)
imagesc(cmp2,image_neptune_clear+image_neptune_stf2);
imagesc(cmp3,image_neptune_overlapped);








