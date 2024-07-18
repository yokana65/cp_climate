% Markus L Fischer
% 2024/03

% Importing KL15 data for TSA
clear, clc, close all

% get subfolders
folder = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(folder));


% load the data
script_read_data_200kyr_all

% display the data
script_display_200ka