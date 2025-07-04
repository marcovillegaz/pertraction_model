function LnGamma = UNIFAC_vdW_FV(compoundsLibrary,unifacData,x,T,n)
% This function calculates the activity coefficient of every component in a
% misxture using the UNIFAC method described in:
% Poling, B. E., Prausnitz, J. M., & O’Connell, J. P. (2001). 
% The Properties of Gases and Liquids (5th ed.). McGraw-Hill.
% Input:
%   data: cell that contains the groups, and the number of groups in
%   compounds. it follows the next structure.
%         data = {
%              group   #1 #2       R     Q    v1   v2 ..... vn
%             'CH3'    1  1   0.9011 0.848    1    1
%             'CH2'    1  2   0.6744 0.540    0    3
%             'CH3CO'  9  18  1.6724 1.488    1    0
%                }
%   aij_Matrix: matrix that contains the group interaction parameters.
%   rho: density of components in [g/cm3]    [rho(1);rho(2);...;rho(n)]
%   MW: moleculas weight of components       [MW(1);MW(2);...,MW(n)]
%   x: molar fraction of both components [x1;x2;...;xn]
%   n: degree of polimerization
%   T: absolute temperatura of the system [k]
% Ouput:
%   Gamma: activity coefficient vector [gamma_1;gamma_2;...;gamma_n]

%% 1) Extracting data
x = x';

% Redefine unifac data
Aij = unifacData.Aij     % group interaction parameters
groupData = unifacData.groups  % cell with groups data

% get component lists
component = fields(compoundsLibrary);

% Extract molar weight and density as array
MW = extractPropertyAsArray(compoundsLibrary, {"molarWeigth",1})
rho_fh = extractPropertyAsFunction(compoundsLibrary, "density")

% MW = data(:,[1 6]);       % molar weight of group
groupData(:,6) = [];             % deleting MW column

% CONSTANTS
[g,c] = size(groupData);  % g = # of groups
c = c-5             % c = # of components


v = MW./rho_fh(T);    % molar volume vector [cm3/mol]

z = 10;     % UNIFAC model constnat
%% 2) R and Q calculations
% Data reaarrengment
r = zeros(c,1);
q = zeros(c,1);
l = zeros(c,1);
aux = cell2mat(groupData(:,4:end));

% r, q and l calculation
for k = 1:c
    r(k) = sum(aux(:,1).*aux(:,k+2));
    q(k) = sum(aux(:,2).*aux(:,k+2));
    l(k) = (z/2)*(r(k)-q(k))-(r(k)-1);
end
    
%% 3) Combinatorial part
    r(end) = n*r(end);
    % normal phi and theta calculation
    phi = zeros(c,1);
    theta = zeros(c,1);
    for k = 1:c
        phi(k) = (r(k)*x(k))/sum(r.*x);
        theta(k) = (q(k)*x(k))/sum(q.*x);
    end
    
    % The activity is calculated for the componentes 'in' the polymer
    LnGamma_C = zeros(c,1);
    for k = 1:c
        A = log(phi(k)/x(k));
        B = (z/2)*q(k)*log(theta(k)/phi(k));
        C = -(phi(k)/x(k))*sum(x.*l);
        LnGamma_C(k) = A + B + C + l(k);
    end

%     % phi, phi_x and theta calculation. (phi_x = phi') by Zhong
%     % Volume fraction correction (Zhong, 1996)
%     r(end) = n*r(end)
%     r_zhong = r
%     r_zhong(end) = r_zhong(end)*0.6593
%     
%     phi = zeros(c,1);
%     phi_x = zeros(c,1);   % (phi')
%     theta = zeros(c,1);
%     for k = 1:c
%         phi(k) = (r(k)*x(k))/sum(r.*x)
%         phi_x(k) = (r_zhong(k)*x(k))/sum(r_zhong.*x)
%         theta(k) = (q(k)*x(k))/sum(q.*x)
%     end
    
%     % The activity is calculated for the componentes 'in' the polymer
%     LnGamma_C = zeros(c,1)
%     for k = 1:c
%         A = log(phi_x(k)/x(k)) + 1 - (phi_x(k)/x(k))
%         B = log(phi(k)/theta(k)) + 1 - (phi(k)/theta(k))
%         LnGamma_C(k) = A - (z*q(k)*0.5)*B
%     end

%% 4) Free volume term
v_h = r.*15.17;    % Hardcore volume cm3/mol
v_fv = v - v_h;    % Free volume cm3/mol

phi_fv = zeros(c,1);
phi_h = zeros(c,1);

for k = 1:c
    phi_fv(k) = (x(k)*v_fv(k))/sum(x.*v_fv);
    phi_h(k) = (x(k)*v_h(k))/sum(x.*v_h);
end

LnGamma_FV = zeros(c,1);
for k = 1:c
    LnGamma_FV(k) = log(phi_fv(k)/phi_h(k)) + ((phi_h(k) - phi_fv(k))/x(k));
end

%% 5) Resdiual part (by groups)
    %% Data reaarrengment.
    % Now, rows are groups and columns are compound
    Q = cell2mat(groupData(:,5));
    V = cell2mat(groupData(:,6:end));
    
    % psi matrix
    psi = exp(-Aij./T);

    %% Residual activity coeff. of group k in a reference pure solution
    % theta_m(k,i)  Ln{gamma}^{i}_{k} 
    % (k) = group  (i) = compound

    % Group molar fraction refering the pure component (i)
    X = zeros(g,c);
    for k = 1:g  % group loop
        for i = 1:c % component loop
            X(k,i) = V(k,i)/sum(V(:,i));
        end
    end
    
    % Area fraction of group (capital theta)
    theta_m = zeros(g,c);
    for k = 1:g % group loop
        for i = 1:c % component loop
            theta_m(k,i) = (Q(k)*X(k,i))/sum(X(:,i).*Q);
        end
    end
    
    % Residual activity coeff. of group k in a reference pure solution
    LnGamma_pure = zeros(g,c);
    for i = 1:c % component loop
        for k = 1:g % group loop
            A = log(sum(theta_m(:,i).*psi(:,k)));
    
            aux = zeros(g,1);
            for m = 1:g % group loop
                aux(m) = theta_m(m,i)*psi(k,m)/sum(theta_m(:,i).*psi(:,m));
            end
            B = sum(aux);
            LnGamma_pure(k,i) = Q(k)*(1-A-B);
        end
    end

    %% Residual activity coeff. of group k in mixtura
    % Group mole fraction refering the group mixture solution
    X = zeros(g,1);
    for k = 1:g
        X(k) = sum((x').*V(k,:))/sum((x').*sum(V,1));
    end
    
    % Area fraction of group
    theta_m = zeros(g,1);
    for k = 1:g % group loop
        theta_m(k) = (Q(k)*X(k))/sum(X.*Q);
    end
    
    % Residual activity coeff. of group k in the mixture
    LnGamma_mix = zeros(g,1);
    for k = 1:g % group loop
        A = log(sum(theta_m.*psi(:,k)));
    
        aux = zeros(g,1);
        for m = 1:g % group loop
            aux(m) = theta_m(m)*psi(k,m)/sum(theta_m.*psi(:,m));
        end
        B = sum(aux);
        LnGamma_mix(k) = Q(k)*(1-A-B);
    end

    %% Final redidual term
    LnGamma_R = zeros(c,1);
    for i = 1:c
        suma = 0;
        for k = 1:g
            suma = suma + V(k,i)*(LnGamma_mix(k)-LnGamma_pure(k,i));
        end
        LnGamma_R(i) = suma;
    end


%% 6) Final activity coefficient
LnGamma = LnGamma_C + LnGamma_R + LnGamma_FV;
LnGamma(isnan(LnGamma)) = 0;   
% compunds with molar fraction equal to zero dosen't have
% activity coefficient

end