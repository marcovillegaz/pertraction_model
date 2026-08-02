classdef UNIFACFVModel < UNIFACModel
    % UNIFACFVModel
    % UNIFAC-FV activity coefficient model: combinatorial + residual +
    % free-volume, for polymer-solvent systems. Delegates to unifacFV.m
    % (formerly UNIFAC_test.m). Requires each compound to carry FVP and
    % density data -- see docs/theory/thermodynamic-models/unifac.md.

    properties (Constant)
        RequiredCompoundProps = ["FVP", "density"]
    end

    methods
        function obj = UNIFACFVModel(compoundsLibObj, thermoLibObj, n)
            arguments
                compoundsLibObj
                thermoLibObj
                n (1,1) double = 1000
            end
            obj@UNIFACModel(compoundsLibObj, thermoLibObj, n);
        end

        function LnGamma = computeActivityCoefficient(obj, temperature, molarFraction)
            LnGamma = unifacFV(obj.compoundsLib, obj.thermoLib, ...
                temperature, molarFraction, obj.polymerizationDegree);
        end
    end
end
