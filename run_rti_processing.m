% Clear workspace and command window
clear all;
close all;
clc;

% Specify your data directory
rti_dir = '/work/imvia/ra7916lu/illumi-net/data/subset/2024_02_22_1_3/images/Face_A/rti_sub_images';

% Set parameters
params = struct();
params.dmd_modes = 45;      % Number of DMD modes
params.hsh_order = 3;       % HSH order
params.save_results = true; % Save the results
params.compute_errors = true; % Compute error metrics

% Get total number of images from LP file
lp_file = dir(fullfile(rti_dir, '*.lp'));
fid = fopen(fullfile(rti_dir, lp_file(1).name));
temp_L = textscan(fid, '%s %f %f %f');
fclose(fid);
total_images = length(temp_L{1})-1;  % -1 for header

% Define indices (adjust these based on your needs)
fit_indices = [55, 22, 76, 44, 72, 15, 42, 40, 9, 85, 11, 103, 78, 28, 79, 5, 62, 56, 39, 35, 16, 66, 34, 7, 43, 68, 69, 27, 19, 84, 25, 73, 49, 13, 24, 3, 17, 38, 8, 81, 6, 67, 36, 91, 83, 54, 50, 70, 46, 100, 61, 101, 97, 41, 58, 48, 90, 57, 75, 32, 98, 59, 63, 88, 37, 29, 93, 1, 52, 21, 2, 23, 87, 95, 74, 86, 82, 20, 60, 71, 14, 92, 51, 102];  % Use 80% for fitting
relight_indices = [30, 65, 64, 53, 45, 94, 104, 47, 10, 0, 18, 31, 89, 96, 77, 4, 80, 33, 12, 26, 99];  % Rest for relighting

fit_indices = fit_indices +1 ;
relight_indices = relight_indices +1 ;

% Run the processing
processRTISelective(rti_dir, fit_indices, relight_indices, params);