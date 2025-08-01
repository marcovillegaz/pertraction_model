function Gamma = thermodynamicsFactors(compoundLibrary,unifacLibrary,temperature,x)
% This function computes the thermodinamic factors matrix from the
% Maxwell-Stefan model, following equation (17) from Kubaczka (2014).
% Input:
%   UNIFAC: struct that contains the UNIFAC data necessary to compute the
%   activity coefficients. with fields .Aij and .groups
%   x: molar fraction vector   [x(1);x(2);...;x(n)]
%   T: temperature of the system [K]
% Output:
%   Gamma: thermodynamic factor matrix (m-1) x (m-1)
% 
% Kubaczka, A. (2014). Prediction of Maxwell-Stefan diffusion coefficients 
% in polymer-multicomponent fluid systems. Journal of Membrane Science,
% 470, 389–398. https://doi.org/10.1016/j.memsci.2014.06.055

%% default
h = 0.001;       % step value for finite difference (COULD BE MODIFIED)
m = length(x);   % # of components
n = 1000;  % degree of polimerization (Zhong,1996) (MODIFY)

%% Finite difference Jacobian
J = zeros(m);    % preallocation


% This loop apply forward finite difference to molar fraction
for i = 1:m  
    x_l = x;
    x_u = x;
    x_u(i) = x_u(i) + h;
    
    % This function generates a column vector of the Ln(gamma) of the components
    LnGamma_l = activityCoefficients('unifac-test', compoundLibrary,unifacLibrary, x_l, temperature, n);
    LnGamma_u = activityCoefficients('unifac-test', compoundLibrary,unifacLibrary, x_u, temperature, n);

    % Forward finite difference
    J(:,i) = (LnGamma_u - LnGamma_l)./(h)  
end


%% Thermodynamic factormatrix calculation
Gamma = zeros(m-1); % preallocation
for i = 1:m-1 
    for k = 1:m-1
        kroneckerDelta = double(i == k); % compute the kronecker delta
        Gamma(i,k) = kroneckerDelta + x(i)*(J(i,k) - J(i,m));
    end
end

end