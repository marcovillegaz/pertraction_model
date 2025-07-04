function Gamma = thermodynamicsFactors(S,UNIFAC_data,x,T)
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

%% Finite difference jacobian
J = zeros(m);    % preallocation

for i = 1:m % this loop apply forward finite difference to molar fraction 
    x_l = x;
    
    x_u = x;
    x_u(i) = x_u(i) + h;
    
    % change unifac function  ###
    % The UNIFACs is difficult to differentiati
    % 
    % 
    % function generates a column vector of the Ln(gamma) of
    % the components

    n = 1;  % degree of polimerization (Zhong,1996) ##
    lngamma_l = UNIFAC_vdW_FV(S,UNIFAC_data,x_l,T,n);
    lngamma_u = UNIFAC_vdW_FV(S,UNIFAC_data,x_u,T,n);

    J(:,i) = (lngamma_u - lngamma_l)./(h)   % forward finite difference
end


%% Thermodynamic factormatrix calclation 
Gamma = zeros(m-1); % preallocation
for i = 1:m-1 
    for k = 1:m-1
        Gamma(i,k) = Kronecker_Delta(i,k) + x(i)*(J(i,k) - J(i,m));
    end
end



function [d] = Kronecker_Delta(j,k)
% This function computes the Kronecker Delta

%% default
if nargin < 2
    error('Too few inputs.  See help KronD')
elseif nargin > 2
    error('Too many inputs.  See help KronD')
end

%% Kronecker_Delta
if j == k
    d = 1;
else
    d = 0;
end

end

end