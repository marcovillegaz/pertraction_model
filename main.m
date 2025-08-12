addpath(genpath('source'))

% Constants
COMPOUNDS_LIST = {"benzene","methylAcetate","polystyrene"};  
COMPOUNDS_FOLDER = "data/test-compounds";
% System values (temporal)
temperature = 300;
molarWeight = [0.1, 0.6, 0.3];

% Loading libraries
compoundsLib = CompoundsLibrary(COMPOUNDS_LIST,COMPOUNDS_FOLDER)
unifacLib = UNIFACLibrary("data/unifac-data","unifac-test.xlsx")

% Loading models 
thermoModel = UNIFACModel(compoundsLib,unifacLib)

% Commpute activity coefficient
Lngamma = thermoModel.computeActivityCoefficient(temperature,molarWeight)



% %% COMPUTE SOMETHING
% fickDiffusivity = computeFickDiffusivity(...
%     compoundLibrary,unifacLibrary, temperature, molar_fraction)