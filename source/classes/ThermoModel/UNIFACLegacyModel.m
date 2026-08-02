classdef UNIFACLegacyModel < UNIFACModel
    % UNIFACLegacyModel
    % Legacy UNIFAC activity coefficient model: combinatorial + residual
    % only, no free-volume correction. Suitable for non-polymer systems.
    % Delegates to unifacLegacy.m.

    properties (Constant)
        RequiredCompoundProps = string.empty  % no extra compound data needed
    end

    methods
        function obj = UNIFACLegacyModel(compoundsLibObj, thermoLibObj)
            % Polymerization degree is not meaningful for legacy UNIFAC;
            % pass 1 through to the (unused) base property.
            obj@UNIFACModel(compoundsLibObj, thermoLibObj, 1);
        end

        function LnGamma = computeActivityCoefficient(obj, temperature, molarFraction)
            LnGamma = unifacLegacy(obj.compoundsLib, obj.thermoLib, ...
                temperature, molarFraction, 1);
        end
    end
end
