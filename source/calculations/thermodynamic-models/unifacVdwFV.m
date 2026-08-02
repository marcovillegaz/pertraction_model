function LnGamma = unifacVdwFV(compoundsLibrary, unifacData, T, x, n)
%UNIFACVDWFV UNIFAC-vdW-FV activity coefficient model (combinatorial +
%   residual + free-volume), for polymer-solvent systems. Ports
%   UNIFAC_vdW_FV.m onto the current CompoundsLibrary/UNIFACLibrary API
%   and the shared unifac/ terms.
%
%   LnGamma = unifacVdwFV(compoundsLibrary, unifacData, T, x, n)
%
%   *** NEEDS REFERENCE ***
%   The original UNIFAC_vdW_FV.m computed its free-volume term with the
%   exact same formula as UNIFAC-FV (unifacFreeVolume.m: v_h = r*15.17,
%   Bondi hard-core volume). "vdW-FV" implies a van der Waals-volume-based
%   free-volume correction distinct from that -- but no such distinction
%   exists anywhere in this repository's source or docs. This function
%   currently reuses unifacFreeVolume.m unchanged, so it is numerically
%   IDENTICAL to unifacFV.m. It is safe to call but is not yet a distinct,
%   validated model. Replace the free-volume term below once the
%   literature reference for UNIFAC-vdW-FV is supplied, and update
%   docs/theory/thermodynamic-models/unifac.md accordingly.
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
%     LnGamma - Ln(activity coefficient) column vector
%
% Reference: Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
%   The Properties of Gases and Liquids (5th ed.). McGraw-Hill, ch. 8.10.
%   [vdW-FV free-volume term: reference pending]

params = unifacGroupParams(compoundsLibrary, unifacData, T, n);

LnGamma_C = unifacCombinatorial(params, x);
LnGamma_R = unifacResidual(params, x);
LnGamma_FV = unifacFreeVolume(params, x);   % TODO: replace with true vdW-FV term

LnGamma = LnGamma_C + LnGamma_R + LnGamma_FV;
LnGamma(isnan(LnGamma)) = 0;

end
