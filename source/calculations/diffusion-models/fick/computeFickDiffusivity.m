function fickDiffusivity = computeFickDiffusivity(...
    compoundLibrary, thermoModel, temperature, molar_fraction)
%MAXWELLSTEFANDIFFUSSION This function compute the effective difussivity,
%that led us work with the mass transfer equation from the Maxwell stefan
%form to the FIcks law form.
%   input:
%       compoundLibrary (CompoundsLibrary): compounds in the system and
%       all their properties.
%       thermoModel (ThermoModel): activity-coefficient model (e.g. a
%       UNIFACFVModel) used to build the thermodynamic factor matrix.
%       The model already owns its unifacLibrary and any model-specific
%       options (e.g. degree of polymerization).
%       temperature: temeprature of the system [K]
%       molar_fraction (array): molar fraction of componentes in the membrane. The
%       length of this array must be the same as the number of compounds in
%       compoundLibrary.
%   Output:
%       DifussionCoefficient (array): n by n matrix that contain all the
%       diffusion coefficientes acoording to the Maxwell-Stefan Diffusion
%       Theory. [cm2/s]
%
% Reference of Maxwel-Stefan theory:
%   Taylor, R., & Krishna, R. (1993). Multicomponent Mass Transfer (1st ed.).
%       John Wiley & Sons.

%% CODE
% Mutual diffusion or Maxwel-Stefan diffusion matrix [-D-]
debugMsg("Mutual difussion calculation")
mutual_diffusion = mutualDiffusion(compoundLibrary,temperature,molar_fraction)
% disp(B)

% Inverse drag coefficient matrix [B]
B = Bmatrix(mutual_diffusion,molar_fraction)
% fprintf('\n%-30s\n\n',"Diffusion coefficient matrix [B] ... OK! ")
% disp(B)

% Thermodynamic factor matrix $[\Gamma]$
Gamma = thermodynamicsFactors(thermoModel,temperature,molar_fraction)
%fprintf('%-30s\n\n',"Thermodinamic factor matrix [Gamma] ... OK! ")
%disp(Gamma)

% Effective diffusion or Fick's diffusion matrix
fickDiffusivity = B\Gamma;   % (cm2/s)

end

