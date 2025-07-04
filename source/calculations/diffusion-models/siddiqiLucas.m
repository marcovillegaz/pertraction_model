function D0_infinite = siddiqiLucas(temperature, compoundLibrary)
% This function calculates the diffusion coefficient at infinite dilution
% for polar and nonpolar systems according to Siddiqi&Luqas correlation. 
% This function only apply for non polymer compounds. 
%
% Inputs:
%   T              - Temperature [K]
%   compoundNames  - Cell array of compound names (strings)
%   boilingVolumes - Array of molar volumes at boiling point [cm3/mol]
%   viscosityFuncs - Cell array of function handles for viscosity [Pa·s]
%
% Output:
%   D0_infinite    - (n)x(n) matrix with diffusion coefficients for the non
%                       polymer interactions [cm2/s]
% 
% Siddiqi, M. A., & Lucas, K. (1986). Correlations for prediction of 
% diffusion in liquids. The Canadian Journal of Chemical Engineering, 
% 64(5), 839–843. https://doi.org/10.1002/cjce.5450640519

%% Extraction of parameters from compoundLibrary
compound_names = fieldnames(compoundLibrary);   % field names

viscosityFuncs = extractPropertyAsArray(compoundLibrary, {"Viscosity"})
% Change this to the real property path
boilingVolumes = extractPropertyAsArray(compoundLibrary, {"boilingPoint",1})

T = temperature;

%% Code
n = length(compound_names);
D0_infinite = zeros(n-1);

for i = 1:n-1 % solute loop
    Vb_1 = boilingVolumes(i); % solute molar volume

    for j = 1:n-1 % solvent loop
        if strcmpi(compound_names{j}, 'water')
            mu_w = viscosityFuncs{j}(T) * 1000; % water viscosity [mPa·s]
            D0_infinite(i,j) = 2.98e-7 * mu_w^(-1.026) * Vb_1^(-0.5473) * T;
        else
            mu_2 = viscosityFuncs{j}(T) * 1000; % solvent viscosity [mPa·s]
            Vb_2 = boilingVolumes(j);           % solvent molar volume
            D0_infinite(i,j) = 9.86e-8 * mu_2^(-0.907) * Vb_1^(-0.45) * Vb_2^(0.265) * T;
        end
    end
end

end