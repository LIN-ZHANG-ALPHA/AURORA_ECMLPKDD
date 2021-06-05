
clear all;  close all; % Clear all previous variables and close all previous windows.

addpath(genpath(pwd))

load('sample_data.mat')

%% Dictionary Parameters

% ++++++++++ periodic Dictionary
Pmax            = 75; %The largest period spanned by the NPDs
Dictionary_pool = {'Ramanujan','NaturalBasis','random' };%Type of the dictionary
Dictionary_type = Dictionary_pool{1};
opts.Dictionary_type = Dictionary_type;
opts.Pmax               = Pmax;

% ++++++++++ spline Dictionary
degree         = 3; % cubic, set as default
num_space = 20; % for equally spaced knot vector
opts.A         = construct_spline_dictionary(1:size(sample_data,1),num_space,degree);

%% 
opts.lambda_1        = .1;
opts.lambda_2        = .1;
opts.lambda_3        = .1;

opts.rho             = 1e-3;
opts.visual          = 0;
opts.max_iter     = 100;

%% Model Parameters

% main function 
[Output_ours,running_time]= AURORA(sample_data,opts);




