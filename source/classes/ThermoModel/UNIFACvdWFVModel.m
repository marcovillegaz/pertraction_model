classdef UNIFACvdWFVModel < UNIFACModel
    % UNIFACvdWFVModel
    % UNIFAC-vdW-FV activity coefficient model: combinatorial + residual +
    % free-volume, for polymer-solvent systems. Delegates to
    % unifacVdwFV.m (formerly UNIFAC_vdW_FV.m).
    %
    % *** NEEDS REFERENCE *** -- see unifacVdwFV.m. The free-volume term
    % currently used here is numerically identical to UNIFACFVModel's;
    % the literature distinction for "vdW-FV" has not yet been supplied.
    % Requires each compound to carry FVP and density data -- see
    % docs/theory/thermodynamic-models/unifac.md.

    properties (Constant)
        RequiredCompoundProps = ["FVP", "density"]
    end

    methods
        function obj = UNIFACvdWFVModel(compoundsLibObj, thermoLibObj, n)
            arguments
                compoundsLibObj
                thermoLibObj
                n (1,1) double = 1000
            end
            obj@UNIFACModel(compoundsLibObj, thermoLibObj, n);
        end

        function LnGamma = computeActivityCoefficient(obj, temperature, molarFraction)
            LnGamma = unifacVdwFV(obj.compoundsLib, obj.thermoLib, ...
                temperature, molarFraction, obj.polymerizationDegree);
        end
    end
end
