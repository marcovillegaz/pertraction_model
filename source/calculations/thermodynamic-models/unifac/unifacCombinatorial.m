function LnGamma_C = unifacCombinatorial(params, x)
%UNIFACCOMBINATORIAL UNIFAC combinatorial (entropic) activity term.
%   LnGamma_C = unifacCombinatorial(params, x), where params is the
%   struct returned by unifacGroupParams.m and x is the molar fraction
%   vector.
%
% Reference: Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
%   The Properties of Gases and Liquids (5th ed.). McGraw-Hill, ch. 8.10.

x = x(:);
c = params.c;
r = params.r;
q = params.q;
l = params.l;
z = params.z;

phi = zeros(c,1);
theta = zeros(c,1);
for k = 1:c
    phi(k) = (r(k)*x(k))/sum(r.*x);
    theta(k) = (q(k)*x(k))/sum(q.*x);
end

LnGamma_C = zeros(c,1);
for k = 1:c
    A = log(phi(k)/x(k));
    B = (z/2)*q(k)*log(theta(k)/phi(k));
    C = -(phi(k)/x(k))*sum(x.*l);
    LnGamma_C(k) = A + B + C + l(k);
end
end
