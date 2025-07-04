
addpath(genpath('source'))

% Define compound to load.
polymerList = {"polystyrene"}
compoundNames = {"benzene","methylAcetate","polystyrene"}  
compoundDataFolder = "data/test-compounds";

x = [0.1, 0.6, 0.3]    % test molar fraction
T = 273.15
n=10000

%% LOAD COMPOUND DATA
% Load compound data and save in library
compoundLibrary = initCompoundLibrary(compoundNames,compoundDataFolder) 
% Load UNIFAC DATA
unifacLibrary = loadUnifacData("unifac-test.xlsx")

%% COMPUTE SOMETHING
% Compute activity coefficient
% LnGamma = activityCoefficients('unifac-vdw-fv', compoundLibrary,unifacLibrary, x, T, n)
LnGamma = activityCoefficients('unifac-test', compoundLibrary,unifacLibrary, x, T, n)


% % Save to no reload every time
% save(fullfile(compoundDataFolder,"compoundLibrary.mat"),"compoundLibrary") 
% 
% load(fullfile(compoundDataFolder,"compoundLibrary.mat"))
% % Compute mutual diffusion ciefficient matrix
mutual_diffusion = mutualDiffusion(x,300,compoundLibrary)


