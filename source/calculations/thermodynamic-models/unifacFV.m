function LnGamma = unifacFV(compoundsLibrary, unifacData, T, x, n)
%UNIFACFV UNIFAC-FV activity coefficient model (combinatorial + residual
%   + free-volume), for polymer-solvent systems. Supersedes UNIFAC_test.m,
%   which was a production model despite the "_test" name.
%
%   LnGamma = unifacFV(compoundsLibrary, unifacData, T, x, n)
%
%   Inputs:
%     compoundsLibrary - CompoundsLibrary object
%     unifacData         - UNIFACLibrary object
%     T                   - absolute temperature [K]
%     x                   - molar fraction vector [x(1);x(2);...;x(n)]
%     n                   - degree of polymerization of the polymer
%                           component (last compound in the library)
%
%   Output:
%     LnGamma - Ln(activity coefficient) column vector [LnGamma_1;...;LnGamma_n]
%
%   Composed from the shared UNIFAC terms in unifac/:
%     unifacGroupParams, unifacCombinatorial, unifacResidual, unifacFreeVolume.
%   See also: unifacLegacy.m, unifacVdwFV.m.
%
% Reference: Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
%   The Properties of Gases and Liquids (5th ed.). McGraw-Hill, ch. 8.10.

params = unifacGroupParams(compoundsLibrary, unifacData, T, n);

LnGamma_C = unifacCombinatorial(params, x);
LnGamma_R = unifacResidual(params, x);
LnGamma_FV = unifacFreeVolume(params, x);

LnGamma = LnGamma_C + LnGamma_R + LnGamma_FV;
LnGamma(isnan(LnGamma)) = 0;
% compounds with molar fraction equal to zero don't have an activity coefficient

end
