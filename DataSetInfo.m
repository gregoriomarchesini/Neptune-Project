% Developer : Gregorio Marchesini
% Date      : 17/02/2021
% Note      : The function 'fitsreader' was developed by professor Lorenz Roth
%             and it is used through out the sript to read fits files


% Description of the observation 
% ––––––––––––––––––––––––––––––

% This script is used to display and manage FITS files obtained from the
% observing campaign of Neptune by the Hubble Space Telescope (HST). The
% observations are made throught the Ultraviolet Camera mounted on HST and
% the target radiation is the Lyman-alpha emission form Neptune.
% The desired emission will be obtained by means of image
% subtratction. In particular the throughput of the UAV filter F25LYA, by means 
% the Lyman-alpha emission could be obtained from direct observation of the
% planet, is very low. A valid altrenative is to subtract a clear image
% (no filter) of the planet and a second image obtained with a F25STF2
% filter which rejects all the Lyman-alpha radiation. The final throughput
% deriving from this subtraction is higher then is the previous approach.


%% Initialization

clc
clear all
close all

%% Target Images from the folder

dir_list   = dir(fullfile(pwd,'HST'));                    % Find all files with the start name 'odq'
target_dir = fullfile({dir_list.folder},{dir_list.name}); % The target files are insife dir_list in a structure.
                                                          % note that each FITS file has the same name of the                                                       
                                                          % directory in which it is contained apart form teh estension           
Observation_name = {dir_list.name}; 
mask             = ~strncmp(Observation_name,'.',1);      % eliminates the first two useless files '..' ans '.' 
Observation_name = Observation_name(mask);
target_dir       = target_dir(mask);                      % The target names of the folder
%% Open the images

fprintf('Loading Files ...\n\n')

for i=1:length(Observation_name)
    [info{i},set_DATA{i}] = fitsreader(target_dir{i},Observation_name{i});
end

%% IMPORTANT : Info from the _apt file. In this files you have calibration informations

folderpath = (fullfile(pwd,'additional_FITS_files'));
dir_int    = dir(folderpath);
names      = {dir_int.name};
index      = find(contains(names,'_apt'));


%% Through Put Informations Stored in the Fits File 

path_useful     = fullfile(dir_int(index).folder,dir_int(index).name);
throughput_info = fitsinfo(path_useful);

% Save the informations on a Text File 

% fileID                 = fopen('Throughputs.txt','w');
% throughput_dataTable   = throughput_info.PrimaryData.Keywords;
% throughput_Binarytable = throughput_info.BinaryTable.Keywords;
% 
% fprintf(fileID,'%s\n',throughput_dataTable{:,3});
% fprintf(fileID,'%s\n',throughput_Binarytable{:,3});
% fclose(fileID);

fprintf('Finished Loading Files ...\n')
close all


% Further Developments :
% ––––––––––––––––––––––

% The script can be converted into a function which will be able to open
% all the files in a given folder just giving the name of the folder.


