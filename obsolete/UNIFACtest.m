function LnGamma = UNIFACtest(UNIFAC_data,x,T)
% This function calculates the activity coefficient of every component in a
% misxture using the UNIFAC method described in:
% Poling, B. E., Prausnitz, J. M., & O’Connell, J. P. (2001). 
% The Properties of Gases and Liquids (5th ed.). McGraw-Hill.
% Input:
%   UNIFAC.groups: cell that contains the groups, and the number of groups in
%   compounds. it follows the next structure.
%         data = {
%              %  group     #1   #2        R      Q        MW   v1   v2
%                 'CH3'      1    1   0.9011  0.848     15.03    0    2
%                 'CH2'      1    2   0.6744  0.540     14.03    0    1
%                 'C'        1    4   0.2195  0.000     12.03    0    1
%                 'ACH'      3    9   0.5313  0.400     13.03    6    0
%                }
%   UNIFAC.Aij: matrix that contains the group interaction parameters.
%   x: molar fraction of both components [x1;x2;...;xn]
%   T: absolute temperatura of the system [k]
% Ouput:
%   Gamma: activity coefficient vector [gamma_1;gamma_2;...;gamma_n]

%% Extracting data
Aij = UNIFAC_data.Aij;      % group interaction parameters
data = UNIFAC_data.groups;  % cell with groups data

% MW = data(:,[1 6]);       % molar weight of group
data(:,6) = []              % deleting MW column

% constants
[g,c] = size(data);  % g = # of groups
c = c-5;             % c = # of components
z = 10;

%% Combinatorial part
    % Data rearrengment
    r = zeros(c,1);
    q = zeros(c,1);
    l = zeros(c,1);

    aux = cell2mat(data(:,4:end));

    % Compounds loop
    for k = 1:c
        r(k) = sum(aux(:,1).*aux(:,k+2));
        q(k) = sum(aux(:,2).*aux(:,k+2));
        l(k) = (z/2)*(r(k)-q(k))-(r(k)-1);
    end

    % psi and theta calculation
    phi = zeros(c,1);
    theta = zeros(c,1);

    for k = 1:c
        phi(k) = (r(k)*x(k))/sum(r.*x);
        theta(k) = (q(k)*x(k))/sum(q.*x);
    end

    % Combinatorial term
    LnGamma_C = zeros(c,1);
    for k = 1:c
        A = log(phi(k)/x(k));
        B = (z/2)*q(k)*log(theta(k)/phi(k));
        C = -(phi(k)/x(k))*sum(x.*l);
        LnGamma_C(k) = A + B + C + l(k);
    end
    
%% Resdiual part (by groups)
    % Data reaarrengment.
    % Now, rows are groups and columns are compound
    Q = cell2mat(data(:,5));
    V = cell2mat(data(:,6:end));    % Group ocurrence
    
    % psi matrix
    psi = exp(-Aij./T);

    %% Residual activity coeff. of group k in a reference pure solution i
    % Ln{gamma}^{i}_{k}
    % Group mole fraction refering the pure component
    X = zeros(g,c);
    for i = 1:c % component loop
        for k = 1:g  % group loop
            X(k,i) = V(k,i)/sum(V(:,i));
        end
    end
    
    % Area fraction of group
    theta_m = zeros(g,c);
    for i = 1:c % component loop
        for k = 1:g % group loop
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

    %% Group residual activity coeff.
    % Group mole fraction refering the gruop mixture solution
    X = zeros(g,1);
    for k = 1:g
        X(k) = sum((x').*V(k,:))/sum(V.*(x'),'all');
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

%% Final activity coefficient

LnGamma = LnGamma_C + LnGamma_R;

% gamma = exp(LnGamma);
% gamma(isnan(gamma)) = 0;

%% Print results
% fprintf('%8s\n','UNIFAC')
% fprintf('%8s %8s\n','LnGamma_C','LnGamma_R')
% fprintf('%8.4f  %8.4f\n',LnGamma_C(1),LnGamma_R(1))
end