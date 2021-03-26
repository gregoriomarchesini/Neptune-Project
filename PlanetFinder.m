%% Planet Finder
%  –––––––––––––
%-------------------------------------------------------------------------%
% The script is able to identify the position of the planet Neptune inside
% the images retrived from the HST observational campaign of Neptune.
% The script derives the position of the center of the planet in terms of
% bits and encirle the planet inside a circle which resembles the correct
% circular shape of the planet.
%-------------------------------------------------------------------------%

close all
clc
clf
clear all

%% Folders processing and Files withdrowal 

files1 = dir('filters');             % add Images Path (filters Specifications folder)
files2 = dir('HST');                 % add Images Path ( HST images folder)

% withdraw images from the files

files_name = {files1(3:end).name};                               % The first name is the Strontium fluoride             
addpath(files1(1).folder);                                       % the second is the clear MAMA sensor    
addpath(files2(1).folder);                                       % the second is the clear MAMA sensor    

%     INSERT THE IMAGE HERE 
% –––––––––––––––––––––––––––––––            % Check DataSetInfo for more
%                                              information
% -------------------------------
% LIST OF FILES
% -------------------------------
% odq405acq_flt.fits
% odq405aeq_flt.fits
% odq406tuq_flt.fits
% odq406twq_flt.fits
% odq407a5q_flt.fits
% odq408raq_flt.fits
% odq408rcq_flt.fits
% -------------------------------

observation   = 'odq406twq_flt.fits';
image_neptune = fitsread(observation,'image',1);                 % image darectly from neptune gallery
                                                                 % images  are contained inside HST folder
                                                                 % and you can decide which file to upload                                                                % 

fitsreader(fullfile(files2(1).folder,observation),observation);                                                               
                                                                 
% Graphics definition

real_image           = figure('Position',[0 0 600 400]);
real_image_sub1      = axes(real_image);hold on
real_image_sub1.XLim = [0 1024];
real_image_sub1.YLim = [0 1024];

real_grid1           = imagesc(image_neptune);
xlabel('pixel')
ylabel('pixel')
title('pure image')


%% Convolution Phase
 
% Kernel definition

kernel        =  fspecial('disk',2.5);   % disk convolutional kernel (5x5)
mask          =  kernel~=0;              % find all the point that are non-zero
kernel(mask)  =  2;                      % put the convlutional matrix with all points at value 2
                                         % this is a guess value

% standard square convolution matrix 5x5
% other matrices can be used in the convolution
% possible modifications to the square matrix
% (in my case every non-zero term
% in 1)


convoluted_neptune      = conv2(image_neptune,kernel);                                 % obtained convlution
convoluted_neptune      = convoluted_neptune(length(kernel):end,length(kernel):end);   % eliminate the edge raws and columns coming from covolution
Mean_value_convolution  = mean(mean(convoluted_neptune));                              % find the mean value in the convoluted image
Variance                = mean(mean((convoluted_neptune-Mean_value_convolution).^2));  % find the VAriance in the convoluted image
Standard                = sqrt(Variance);                                              % find the standard deviation

% I make a first order assumption on the possible way to eliminate the
% noise from the figure. I take the mean and I add the standard deviation
% as maximum trash hold of the noise. This eliminates the most of the noise

convoluted_neptune(convoluted_neptune<(mean(mean(convoluted_neptune))+Standard))=0; 

% trim of the edges that are the part of the 
% convolution in which the mask is not completely 
% overlapped with the original
% signal. remeber a convluted signal has dimension
% n+m-1 where n is the first matrix (square) dimenion
% and m is the kernel dimension.
% chek that the size convluted matrix and it must
% be equal to the orginal image after the reduction

% figure definition

conv_image = figure('Position',[0 200 1200 500]);  
axes_conv1 = subplot(121);hold on,legend()
axes_conv2 = subplot(122);hold on

imagesc(axes_conv2,convoluted_neptune)
xlim([0 1024])
ylim([0 1024])
colorbar;



% Now what I can do it is using a bigger kernel in order to find the center of
% the planet. This is a possible idea that can be used
% take a kernel rather big

%% Analysis on the Position of the Planet
% First order approximation
% we don't have an estimate of the background noise but we should work on
% that (look in the documenation of the stis files)


% Integration along the rows and the columns

x_coordinate=sum(convoluted_neptune);    % sum along the x-axis
y_coordinate=sum(convoluted_neptune,2);  % sum along the y-axis

% Bad Edges trimming
% The behaviour of the image at the edge can be seriously noisy and 
% mislead the center location

trim_pix=50; %pixel to be set at zero in the count coordinate counts
mask=(1:length(x_coordinate))<trim_pix | (1:length(x_coordinate))>(length(x_coordinate)-trim_pix);

trimmed_x=x_coordinate;
trimmed_y=y_coordinate;
trimmed_y(mask)=0;
trimmed_x(mask)=0;

% This passage is only for smoothen the data 
smoothing_factor=80;                                 % interesting to check how the smoothing 
                                                     % factor is able to impact the 
                                                     % centre location
           
trimmed_x=movmean(trimmed_x,smoothing_factor); % smooth the data with a movemean
trimmed_y=movmean(trimmed_y,smoothing_factor); % smooth the data with a movemean

% I associate the maximum values with the center location of the planet

[max_intensity_x,x_centre]=max(trimmed_x); % find the maximum
[max_intensity_y,y_centre]=max(trimmed_y); % find the maximum
% the x_centre is the positiona at which the colum integration is maximum
% (the rational is that the column is crossing the whole y_axis of a
% circular star and so it is encapsules the longest segment of the star the
% is a diameter)
% same logic for the y-axis


% plot the results

figure('Position',[500 200 500 500])

% trimmed plot
subplot(121)
hold on
grid on

plot(trimmed_x(),'r','linewidth',4,'DisplayName','x_{coordinate}');
plot(trimmed_y(),'k','linewidth',4,'DisplayName','y_{coordinate}');

plot(x_centre,max_intensity_x,'go','linewidth',8);
plot(y_centre,max_intensity_y,'go','linewidth',8);
plot(axes_conv2,x_centre,y_centre,'go','linewidth',8);
xlabel('pixels')
ylabel('integrated pixel values')
title('column-wise and row-wise integration (trimmed image)')
legend

% not trimmed plot
subplot(122)
hold on
plot(x_coordinate,'r','linewidth',4,'DisplayName','x_{coordinate}');
plot(y_coordinate,'k','linewidth',4,'DisplayName','y_{coordinate}');

xlabel('Raw / Columns')
ylabel('integrated counts values')
title('column-wise and row-wise integration (not trimmed)')
legend



%% Symmetry lines of the planet
% The following plot is not really what I was looking for.
% I hoped the curves where going to bulge a bit more once crossing the
% planet, but it is not the case


symmetry_line_x = convoluted_neptune(x_centre,:);
symmetry_line_y = convoluted_neptune(:,y_centre);

figure()
hold on

plot(cumsum(symmetry_line_y),'b-','linewidth',3,'DisplayName','y_{axis}'); % Integrated counts along the axes of symmetry
plot(cumsum(symmetry_line_x),'k-','linewidth',3,'DisplayName','x_{axis}');

xlabel('pixels')
ylabel('pixel values along the centre axis')
title('integrated pixel counts along the axisline')
legend


%% Disk Averaged Mean Approach
% This section is a bit critical. We are trying to locate the exact
% position of a circe encnpassing the planet. In this case we have an
% estimate of  the planet center given by the previous passages, but we
% miss an approximation of the radious. My first attempt is to find out
% what is the integrated value of the pixels in a gaussin kernel that is
% centered arounf the estimated center. I start with a small kernel and I
% then start to expand the kernel littel by little and I save the value of
% the integrated product of the image and the kernel at that point.
% My point is that the resulting function (integrated value vs kernel
% dimension) must have a minimum difference (derivative) at the location of the radius
% of the planet. This is because, further enlargment of the kernel will
% increase less and less significantly the total integral once the
% only backgroun (outside the palanet rim) is estimated
    
%  %----------------------------------------------------------------------------------%



   % varying radious
   
   rad             = 31:0.5:70.5;          % FIRST FILTER TRIAL
   mean_disk_value = zeros(1,length(rad)); % preallocation
   
   % Graphical objects preallocations
   rectangle          = rectangle(axes_conv2);
   disk_move          = plot(axes_conv1,0,0,'-r','linewidth',3,'DisplayName','pixel count sum');
   differential_count = plot(axes_conv1,0,0,'-b','linewidth',3,'DisplayName','pixel count difference');
   legend()
   
   
   
   
   
   axes_conv1.XLim           = [rad(1),rad(end)];
   axes_conv1.XGrid          = 'on';
   axes_conv1.YGrid          = 'on';
   axes_conv1.XLabel.String  = {'radious'};
   axes_conv1.YLabel.String  = {'one step difference'};
   axes_conv1.Title.String  = {'Disk sum difference'};
   
   %% The Algorith logic
   
   % Now that the center of the planet is found we must find the radious of
   % the planet. In order to do that a cicular naular kernel is created
   % with ever increasing radious. What the programm does is multipling the
   % circular anulus kernel times the counts on the portion of the image
   % that is covered by the circular anulus kernel. The circular kernel is
   % centered at the center position of the planet. The process is repeated
   % for verious radii. The correct radious of the planet is found when the
   % difference between the convlution at one step minus the value of the
   % convoulution at the previous step is maximum in absolute value. 
   % The logic is the following : when I am inside the planet, the
   % convolution will give almost the same result if the variation of the
   % circular kernel dimesion is reasonably small. The significant
   % difference wiill be at the moment the circular kernel is aoutside the
   % planet and the number of counts get close to zero, while the
   % preceeding convlution was still inside the planet.
   
   
   
 % ONLY IF YOU NEED THE GIFT
       
%   filename   = 'Planet_gif.gif';  % create a new file (attention to not overwrite)    
%   del        = 0.1;               % Time between animations frame 
   
   
   for i=1:length(rad)
       
       % FIRST FILTER TRIAL DISK
       
       % 1.5 2.5 3.5  are all radious of the filter
       % so the matrix is 3x3 5x5 7x7

       disk_filter        = circular_kernel(rad(i),rad(i)-2,'rim',200,'stuff',-400); % the circular kernel is weighted so to reduce it to a sort of Wavelet function
       l                  = floor(rad(i));            % this term serves to center the portion of the image which will be convoluted 
       Area               = pi*rad(i)^2;              % area inside the anlus
       mean_disk_value(i) = sum(sum(convoluted_neptune(x_centre-l:x_centre+l,y_centre-l:y_centre+l)...
           .* disk_filter));
       
       % Figure updating on the convolution image
   
       
       if i>2
           differential_count.YData=abs(diff(mean_disk_value(1:i)));
           differential_count.XData=rad(1:i-1);                  %minus one because diff reduces the size of the vector
           axes_conv1.YLim=[-1,max(differential_count.YData)*1.10];
       end
       %grafical objects
       rectangle.Position=[x_centre-l y_centre-l 2*rad(i) 2*rad(i)];
       rectangle.LineWidth=3;
   
       % once you fix the radious of the filter the matrix is of odd dimesions
       % we need to deine the number of elemnts that are at the left and right
       % of the center pixel in this matrix. This is the length of the matrix
       % divided by two and then apply the floor. In our case this is already
       % the direct floor of i  
       
       
       
          
% ONLY IF YOU NEED THE GIFT
      
%       drawnow 
%       frame      = getframe(conv_image);
%       im         = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);   % This creates a matrix of index and color map (cm)
%       
%       if i == 1;
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
%       else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
%       end
%        
%        
       
    
   end

[~,max_radious_index]=max(abs(diff(mean_disk_value)));                      % find the correct radious now
circle=viscircles(axes_conv2,[x_centre,y_centre],rad(max_radious_index));   % plot a circle in the convoluted neptune image
plot(real_image_sub1,x_centre,y_centre,'or','linewidth',10)                 % draw center in the image

fprintf('Outputs:\n')
fprintf('Planet x coordinate :  %g pixel \n',x_centre);
fprintf('Planet y coordinate :  %g pixel \n',y_centre);
fprintf('Estimated radious   :  %g pixels \n',rad(max_radious_index))

