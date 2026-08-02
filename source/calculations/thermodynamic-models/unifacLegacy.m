function LnGamma = unifacLegacy(compoundsLibrary, unifacData, T, x, n)
%UNIFACLEGACY Legacy UNIFAC activity coefficient model (combinatorial +
%   residual only, no free-volume correction). Suitable for non-polymer
%   systems. Ports the standalone obsolete/UNIFACtest.m onto the current
%   CompoundsLibrary/UNIFACLibrary API and the shared unifac/ terms.
%
%   LnGamma = unifacLegacy(compoundsLibrary, unifacData, T, x, n)
%
%   n (degree of polymerization) is accepted only for signature parity
%   with unifacFV.m / unifacVdwFV.m -- it is UNUSED here. Legacy UNIFAC
%   has no polymer volume correction and no free-volume term.
%
%   Inputs:
%     compoundsLibrary - CompoundsLibrary object
%     unifacData         - UNIFACLibrary object
%     T                   - absolute temperature [K]
%     x                   - molar fraction vector [x(1);x(2);...;x(n)]
%     n                   - unused; kept for API parity across variants
%
%   Output:
%     LnGamma - Ln(activity coefficient) column vector
%
% Reference: Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
%   The Properties of Gases and Liquids (5th ed.). McGraw-Hill, ch. 8.10.

arguments
    compoundsLibrary
    unifacData
    T (1,1) double
    x
    n (1,1) double = 1  % unused; API parity with unifacFV.m / unifacVdwFV.m
end

% No polymer volume correction: call unifacGroupParams without n so the
% default (n = 1, a no-op on r(end)) applies.
params = unifacGroupParams(compoundsLibrary, unifacData, T);

LnGamma_C = unifacCombinatorial(params, x);
LnGamma_R = unifacResidual(params, x);

LnGamma = LnGamma_C + LnGamma_R;
LnGamma(isnan(LnGamma)) = 0;

end
