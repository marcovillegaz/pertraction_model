function D0_infinite = SiddiqiLucas(T,S)
% This function calculate the difussion coefficient at infinite dilution
% for polar and nonpolar systems.
% Input: 
%   T: system temperature [K]
%   component: struct that contains the compounds as fields and the group
%   of propeties. Each compounds has the properties gruops as fields.
% Output:
%   D0_infinite: (n-1)x(n-1) matrix that containsd the diffusion 
% coefficients at infinite dilution [cm2/s] for the nonpolymeric compounds
% 
% Siddiqi, M. A., & Lucas, K. (1986). Correlations for prediction of 
% diffusion in liquids. The Canadian Journal of Chemical Engineering, 
% 64(5), 839–843. https://doi.org/10.1002/cjce.5450640519

%% CODE
component = fieldnames(S);
n = length(component);

% Infinite dilution duffusion coefifient
D0_infinite = zeros(n-1);  % [cm2/s]
for i = 1:n-1      % solute loop
    Vb_1 = S.(component{i}).boilingPoint{2,4};   % molar volume at Tb [cm3/mol] (solute)

    for j = 1:n-1  % solvent loop
        if strcmp(component{j},'water')
            mu_w = S.water.viscosity(T)*1000;    % viscosity of water [mPa s] (solvent)
            D0_infinite(i,j) = 2.98*10^(-7)*(mu_w^(-1.026))*(Vb_1^(-0.5473))*T;

        else
            mu_2 = S.(component{i}).viscosity(T)*1000; % viscosity [mPa s] (solvent)
            Vb_2 = S.(component{j}).boilingPoint{2,4}; % molar volume at Tb [cm3/mol] (solvent)
            D0_infinite(i,j) = 9.86*10^(-8)*(mu_2^(-0.907))*(Vb_1^(-0.45))*(Vb_2^(0.265))*T;
            
        end
    end
end

end
