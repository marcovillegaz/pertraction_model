function B = Bmatrix(Dms,x)
% This function computes the B matrix from the maxwell-stefan equation
% system.
% Input:
%   D_MS: Matrix that contains the mutual-difussion coefficient
%   x: Vector that contains the molar fraccion of the components
% Output:
%   B: Matrix that contains the B coefficents of the MS model.
%
% References:
%   Kubaczka, A., Kamiński, W., & Marszałek, J. (2018). Predicting mass 
%       fluxes in the pervaporation process using Maxwell-Stefan diffusion 
%       coefficients. Journal of Membrane Science, 546(July 2017), 111–119.
%       https://doi.org/10.1016/j.memsci.2017.08.074

%% CODE
% Because the general maxweel-Stefan equations has m-1 equations
% The last components(m) is the polymer
m = length(x) ;            % Number of components  
B = zeros(m-1,m-1);

% B_{ii}
for i = 1:m-1
    sum = 0;
        for k = 1:m-1
            if k ~= i
                sum = sum + (x(k)/Dms(i,k));
            end
        end
    B(i,i) = x(i)/Dms(i,end) + sum;
end

% B_{ik}
for i = 1:m-1
    for k = 1:m-1
        if k ~= i
            B(i,k) = -x(i)*((1/Dms(i,k))-(1/Dms(i,m)));
        end
    end
end

end 