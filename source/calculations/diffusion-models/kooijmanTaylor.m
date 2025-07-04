function D_mutual = kooijmanTaylor(x,D0_infinite)
% This function calculates the mutual diffusion coefficient of the none
% polymer compundas assuming that the mutual interactio od the two
% diffusing components are of similar nature as in liquid.
% Input:
%   x: molar fraction array [x(1),x(2),...,x(n)]
%   D_infinite: (n-1)x(n-1) matrix that contains the diffusion 
% coefficients at infinite dilution [cm2/s]
% Output:
%   D_mutual: nxn matrix that contains the mutual diffusion cofficients for
%   the Maxwell-Stefan diffusion model.
%
% Kooijman, H. A., & Taylor, R. (1991). Estimation of Diffusion Coefficients 
% in Multicomponent Liquid Systems. Industrial and Engineering Chemistry 
% Research, 30(6), 1217–1222. https://doi.org/10.1021/ie00054a023

%% CODE
n = length(x);

D_mutual = zeros(n);
for i = 1:n-1   % compound i
    for j = 1:n-1 % compound j
        sum = 0;
        for k = 1:n-1
            if k ~= i && k ~= j
                aux1 = (x(k)/2)*(log(D0_infinite(i,k))+log(D0_infinite(j,k)));
                sum = sum + aux1;
            end
        end
        aux2 = x(i)*log(D0_infinite(i,j))+x(j)*log(D0_infinite(j,i))+sum;
        D_mutual(i,j) = exp(aux2);
    end 
end

end 