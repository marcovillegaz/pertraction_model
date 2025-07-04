function self_diffusion = vrentasVrentas(molar_fraction, temperature, compoundLibrary)
% This function compute the self difussion coefficients for a multicomponent
% system employin the method of Vrentas and Vrantas derived from the
% free-volume theory. 
% Input:
%   S: struct that contains the compounds as fields and the group
%   of propeties. Each compounds has the properties gruops as fields.
%   x: molar fraction vector [x(1) x(2) ... x(n)]
%   T: system temperature
% Output:
%   D: self difussion vector [D(1),D(2),...,D(n)]
%
% References:
%   Vrentas, J. S., & Vrentas, C. M. (1998). Predictive methods for 
%       self-diffusion and mutual diffusion coefficients in polymer-solvent 
%       systems. European Polymer Journal, 34(5–6), 
%       797–803. https://doi.org/10.1016/S0014-3057(97)00205-X
%   Kubaczka, A., Kamiński, W., & Marszałek, J. (2018). Predicting mass 
%       fluxes in the pervaporation process using Maxwell-Stefan diffusion 
%       coefficients. Journal of Membrane Science, 546(July 2017), 111–119.
%       https://doi.org/10.1016/j.memsci.2017.08.074

%% CODE
% constant
molar_fraction = molar_fraction';
R = 1.987   ;   % gas constant [cal/mol K]

component = fieldnames(compoundLibrary);   % field names
n = length(component);        % number of components

T = temperature; % Redefinition of temperature for visualization in formula
% %% Extraction of molar weigth
% MW = zeros(1,n); % molar weight   [g/mol]
% for i = 1:n
%     MW(i) = S.(component{i}).molarWeigth{1,4};  % g/mol
% end

%% Extraction of parameters from compoundLibrary
% --- Extract molar weight ---
MW = extractPropertyAsArray(compoundLibrary, {"molarWeigth",1}) 

% --- FREE VOLUME PARAMTERS EXTRACTION ---
V_h =  extractPropertyAsArray(compoundLibrary, {"FVP",2})   % molar critical hole [cm3/mol]
D0 =  extractPropertyAsArray(compoundLibrary, {"FVP",3})    % D0_i [cm2/s] 
beta =  extractPropertyAsArray(compoundLibrary, {"FVP",4})  % K1j/gamma [cm3/g K]
delta =  extractPropertyAsArray(compoundLibrary, {"FVP",5}) % K_2i - T_gi [K]
% xi =  extractPropertyAsArray(compoundLibrary, {"FVP",6})    % xi ratio [-]
E =  extractPropertyAsArray(compoundLibrary, {"FVP",7})     % activaiton energy [cal/mol]

xi = V_h./V_h(end) % Compute the jumping unit ratio (xi) [-]

%% Computation
% molar fraction to weight fraction
w = (molar_fraction.*MW)./sum(molar_fraction.*MW);   

% IN FUTURE EXTEND TO WELL KNOWN POLYMERS (Vrenta & Vrentas 1997)
% V_FH/gamma parameter (denominador) According to classic notation
V_FH_gamma = 0;  % [cm3/g]
for i = 1:n
    aux =  w(i).*beta(i).*(delta(i) + T);
    V_FH_gamma = V_FH_gamma + aux;
end

% self difusion calculation for the non polymer component
self_diffusion = zeros(1,n);
for i = 1:n
    suma0 = 0;
    for j = 1:n  % numerador 
        suma0 = suma0 + w(j)*V_h(j)*(xi(i)/xi(j));  % [cm3/g]
    end
    self_diffusion(i) = D0(i)*exp(-E(i)/R*T)*exp(-suma0/V_FH_gamma) % [cm2/s]

end

end