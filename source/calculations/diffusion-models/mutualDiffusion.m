function mutual_diffussion = mutualDiffusion(...
    molar_fraction, temperature, compoundLibrary)
%MUTUALDIFFUSION Summary of this function goes here
%   This function computes the mutual diffusion matrix for n compounds. 
% In the case of mass transfer involving membranes, the last row and column
% correspond to polymer interaction. If other mass transfer process are
% being modeled with the Maxwell-Stefan approuch this part of the code
% should be modified. 
%
%   input:
%       molar_fraction: molar fraction array of compounds
%       T: temperature [K]
%       compoundLibrary_S (struct): that contains all the information of
%       compounds
%       systemCondition_S (struct): contains system information. 

% In future maybe add a condition to differentite between polymer and non
% polymer interaction, this could expand the applicability of this code,
% for example in composite membranes. 

%% Non polimer iteractions
% Compute diffusion coefficients at infinite dilution 
fprintf('\n%-30s\n',"infinite dilution diffusion coefficients")
inf_dilu_diff = siddiqiLucas(temperature, compoundLibrary)

% Compute mutual diffusion coefficients for non polymer interactions
fprintf('\n%-30s\n',"mutual diffusion matrix for non polymer interaction")
mutual_diff_npi = kooijmanTaylor(molar_fraction,inf_dilu_diff)

%% Polymer interaction
% Self diffusion coefficient with Free Volume Theory
fprintf('\n%-30s\n',"Self diffusion coefficients")
self_diff = vrentasVrentas(molar_fraction,temperature,compoundLibrary) 

% Mutual diffusion coefficient for polimer interactions}
fprintf('\n%-30s\n',"Mutual diffusion coefficient for polymers pair")
mutual_diff_pi = kubaczka(molar_fraction,self_diff,mutual_diff_npi)

%% Combine both interaction into one matrix
% Both mutual diffussion coefficient (for polymer and non polymer
% interaction) return as a matrix that can be sum to obtain the mutual
% difussion coefficient matrix. 

mutual_diffussion = mutual_diff_pi
end
