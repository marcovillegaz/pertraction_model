function LnGamma_FV = unifacFreeVolume(params, x)
%UNIFACFREEVOLUME UNIFAC-FV free-volume correction term, used for
%   polymer-solvent systems (Oishi & Prausnitz type free-volume term).
%   LnGamma_FV = unifacFreeVolume(params, x), where params is the struct
%   returned by unifacGroupParams.m (params.r already includes any
%   polymer volume correction, and params.v is the molar volume derived
%   from density) and x is the molar fraction vector.
%
% Reference: Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
%   The Properties of Gases and Liquids (5th ed.). McGraw-Hill, ch. 8.10.

x = x(:);
c = params.c;
r = params.r;
v = params.v;

v_h = r.*15.17;       % Hard-core volume [cm3/mol]
v_fv = v - v_h;        % Free volume [cm3/mol]

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
end
