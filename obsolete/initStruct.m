function compoundLibrary = initStruct(compoundsfiles,UNIFACfile)
% INITSTRUCT Initializes the full system structure from Excel
% DEscription: This function import the properties requiered to modelate the mass
% trasnfer in a membrane using Maxwell-Stefan model for diffusion.
%
% Input:
%   compoundsFIles (cell): cell with the filename of  each compounds.
%   UNIFACFile (str): Name of the UNIFAC file.
% Output:
%   compoundLibrary (struct): object that contain all compounds
%   infromation. This master struct feed the mass trasnfer model. 
 
%
% The experimental data of density and viscosity is fitted and the
% correspondin equations is stored in the corresponding field in the main
% struct. 
%
% Group contribution methos are used for compounds without experimental 
% data. Sastri & Rao for viscosity and GCVOL60 (elbro) for density.

%% Initialize
addpath(genpath('source'));  % Add all subfolders in source

%% IMPORTING AND PREPROCESSING DATA BLOCK
% Importing porperties matrix
filename = 'PCB77_perstract.xlsx';
sheets = {'Properties','UNIFAC'};
ranges = {'A1:G21','A3:O8'};

[S,propGroups,UNIFAC_data] = propertiesImporter(filename,sheets,ranges);

%% Tyn & Calus molar volume at boiing point
S = TynCalus(S)

%% GCVOL density estimation
filename = 'GCVOL.xlsx';
range = ('A1:E7');
S = GCVOL60(S,filename,range);

%% Sastri & Rao viscosity estimation
filename = 'SastriRao.xlsx';
range = ('A1:D6');
S = sastriRao(S,filename,range);

%% Clear and save
clear filename range ranges sheets T_range plotNames 
save compoundData.mat S propGroups UNIFAC_data

end












