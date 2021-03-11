% Developer : Gregorio Marchesini
% Date      : 17/02/2021
% Note      :  The function 'fitsreader' was developed by professor Lorenz Roth
%              and it is used through out the sript to read fits files


% Description of the observation : 
% This script is used to plot and manage FITS files obtained from the
% observing campaign of Neptune by the Hubble Space Telescope (HST). The
% observations are made throught the Ultraviolet Camera mounted on HST and
% the target radiation is the Lyman-alpha emission form Neptune.
% The given desired emission will be obtained by means of image
% subtratction. In each couple of images, one image is obtained with no
% filters in the ultraviolet wavelenght, while the second is obtained
% using a filter that rejects the particular wavelength of the Lyman-alpha
% emission (121 nm). Once subtracting the images it will be possible to
% appricciate only the required wavlength. This technique seems to achieve
% higher results then direct filtration of the spectrum to obtain a direct
% image of the Lyman-alpha emission form Neptune.


%% Initialization

clc
clear all
close all

%% Target Images from the folder
target_folder=fullfile('.','HST');                     % target directory 
dir_list=dir([target_folder,'/','*odq*']);             % Find all files with the start name 'odq'
target_dir=fullfile({dir_list.folder},{dir_list.name}); % The target files are insife dir_list in a structure.
                                                       % note that each FITS file has the same name of the                                                       
                                                       % directory in which it is contained apart form teh estension                                                    % Using {} puts each field of the structure in a cell-arry
Observation_name={dir_list.name};                      % The target names of the folder
%% Open the images
for i=1:length(target_dir)
    [info{i},set_DATA{i}]=fitsreader(target_dir{i},Observation_name{i});
end

%% IMPORTant Manoeuvre to Refind the Throughput Info from the _apt file

cd referencefiles  
dir_int=dir('*_apt*.fits');
cd ..
path_useful=fullfile(dir_int(1).folder,dir_int(1).name);
throughput_info=fitsinfo(path_useful);

fileID=fopen('Throughputs.txt','w');
throughput_dataTable=throughput_info.PrimaryData.Keywords;
throughput_Binarytable=throughput_info.BinaryTable.Keywords;
fprintf(fileID,'%s\n',throughput_dataTable{:,3});
fprintf(fileID,'%s\n',throughput_Binarytable{:,3});
fclose(fileID);


close all




