%% Script for simulation and reconstruction of oligomers

% author:  Magdalena Schneider
% date:    08.08.2020
% version: 1.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add folder and subfolders to path
folder = fileparts(which(mfilename));
addpath(genpath(folder));

% Set random seed
rng_seed = rng("default");

% Saving options
doSaveResults = true;
filename = 'simulatedOligomers';
pathSave = 'testdata/';
if ~exist(pathSave, 'dir')
    mkdir(pathSave)
end


%% Set simulation parameters

% Oligomers
params.oligoDegree = 4;
params.sizeSpecification = 'sidelength';
params.sidelength = 5;
params.numberOligomers = 500000;

% Labeling
params.labelEff = 1;

% Brightness
params.maxIntensity = 100000; % N_max
params.locPrecThreshold = 10;

% Background
params.background = 0; % background noise
params.pxSize = 100;   % pixel size
params.sigma = 160;    % width of PSF

% Blinking statistics
load('blinkDist_logNormal_meanlog5_stdlog2.mat')
params.blink_dist{1}=blink_dist;

% Polarization vectors
pol_azimuth = [0; pi/2]; % orthogonal to each other
pol_elevation = [0; 0];  % orthogonal to optical axis
distance = [1; 1];       % normalized
[excit_x,excit_y,excit_z] = sph2cart(pol_azimuth,pol_elevation,distance);
params.pol_vecs_cart = [excit_x,excit_y,excit_z];


%% Run simulation

disp('Running simulation...')
sim_result = simulateOligomers_cryoSMLM(params);

if doSaveResults
    % mat-file
    save([pathSave,filename,'.mat'],'sim_result')
end

