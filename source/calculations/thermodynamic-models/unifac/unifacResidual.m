function LnGamma_R = unifacResidual(params, x)
%UNIFACRESIDUAL UNIFAC residual (group-interaction) activity term.
%   LnGamma_R = unifacResidual(params, x), where params is the struct
%   returned by unifacGroupParams.m and x is the molar fraction vector.
%   Identical for every UNIFAC variant in this codebase (legacy, FV,
%   vdW-FV) -- only the combinatorial and free-volume terms differ.
%
% Reference: Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
%   The Properties of Gases and Liquids (5th ed.). McGraw-Hill, ch. 8.10.

x = x(:);
g = params.g;
c = params.c;
Q = params.Q;
V = params.V;
psi = params.psi;

%% Residual activity coeff. of group k in a reference pure solution i
% Group mole fraction referring to the pure component (i)
X = zeros(g,c);
for k = 1:g
    for i = 1:c
        X(k,i) = V(k,i)/sum(V(:,i));
    end
end

% Area fraction of group in pure component (i)
theta_m = zeros(g,c);
for k = 1:g
    for i = 1:c
        theta_m(k,i) = (Q(k)*X(k,i))/sum(X(:,i).*Q);
    end
end

LnGamma_pure = zeros(g,c);
for i = 1:c
    for k = 1:g
        A = log(sum(theta_m(:,i).*psi(:,k)));

        aux = zeros(g,1);
        for m = 1:g
            aux(m) = theta_m(m,i)*psi(k,m)/sum(theta_m(:,i).*psi(:,m));
        end
        B = sum(aux);
        LnGamma_pure(k,i) = Q(k)*(1-A-B);
    end
end

%% Residual activity coeff. of group k in the mixture
% Group mole fraction referring to the group mixture solution
X = zeros(g,1);
for k = 1:g
    X(k) = sum((x').*V(k,:))/sum((x').*sum(V,1));
end

% Area fraction of group
theta_m = zeros(g,1);
for k = 1:g
    theta_m(k) = (Q(k)*X(k))/sum(X.*Q);
end

LnGamma_mix = zeros(g,1);
for k = 1:g
    A = log(sum(theta_m.*psi(:,k)));

    aux = zeros(g,1);
    for m = 1:g
        aux(m) = theta_m(m)*psi(k,m)/sum(theta_m.*psi(:,m));
    end
    B = sum(aux);
    LnGamma_mix(k) = Q(k)*(1-A-B);
end

%% Final residual term
LnGamma_R = zeros(c,1);
for i = 1:c
    suma = 0;
    for k = 1:g
        suma = suma + V(k,i)*(LnGamma_mix(k)-LnGamma_pure(k,i));
    end
    LnGamma_R(i) = suma;
end
end
