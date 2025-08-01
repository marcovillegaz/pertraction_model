addpath(genpath('source'))

% Define compound to load.
polymerList = {"polystyrene"}
compoundNames = {"benzene","methylAcetate","polystyrene"}  
compoundDataFolder = "data/test-compounds";

molar_fraction = [0.1, 0.6, 0.3]    % test molar fraction in the membrane
temperature = 273.15 % temperature of the system

%% LOAD COMPOUND DATA
% Load compound data and save in library
compoundLibrary = initCompoundLibrary(compoundNames,compoundDataFolder) 
% Load UNIFAC DATA
unifacLibrary = loadUnifacData("unifac-test.xlsx")

%% COMPUTE SOMETHING
fickDiffusivity = computeFickDiffusivity(...
    compoundLibrary,unifacLibrary, temperature, molar_fraction)