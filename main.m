addpath(genpath('source'))

% Constants
COMPOUNDS_LIST = {"benzene","methylAcetate","polystyrene"};
COMPOUNDS_FOLDER = projectPath("data","input","compounds");
% System values (temporal)
temperature = 300;
molarFraction = [0.1, 0.6, 0.3];

% Loading libraries
compoundsLib = CompoundsLibrary(COMPOUNDS_LIST,COMPOUNDS_FOLDER)
unifacLib = UNIFACLibrary(projectPath("data","input","unifac"),"unifac-test.xlsx")

% Loading models -- pick a UNIFAC variant: "unifac" (legacy), "unifac-fv",
% or "unifac-vdw-fv". See docs/theory/thermodynamic-models/unifac.md.
thermoModel = ThermoModel.create("unifac-fv", compoundsLib, unifacLib)

% Commpute activity coefficient
Lngamma = thermoModel.computeActivityCoefficient(temperature,molarFraction)


% %% COMPUTE SOMETHING
% fickDiffusivity = computeFickDiffusivity(...
%     compoundsLib, thermoModel, temperature, molarFraction)