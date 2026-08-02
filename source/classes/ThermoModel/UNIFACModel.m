classdef UNIFACModel < ThermoModel
    % UNIFACModel
    % Abstract base for the UNIFAC family of activity coefficient models
    % (legacy, FV, vdW-FV). Holds the degree of polymerization shared by
    % the polymer-aware variants (ignored by the legacy variant), and
    % delegates the actual calculation to the composed functions in
    % source/calculations/thermodynamic-models/ (unifacLegacy.m,
    % unifacFV.m, unifacVdwFV.m), which are themselves built from the
    % shared terms in .../unifac/.
    %
    % computeActivityCoefficient and RequiredCompoundProps remain
    % abstract here (inherited from ThermoModel) -- only the concrete
    % leaf classes implement them.
    %
    % See also: UNIFACLegacyModel, UNIFACFVModel, UNIFACvdWFVModel,
    % ThermoModel.create.

    properties (Access = public)
        polymerizationDegree = 1000  % degree of polymerization n (unused by legacy)
    end

    methods
        function obj = UNIFACModel(compoundsLibObj, thermoLibObj, n)
            arguments
                compoundsLibObj
                thermoLibObj
                n (1,1) double = 1000
            end
            obj@ThermoModel(compoundsLibObj, thermoLibObj);
            obj.polymerizationDegree = n;
        end
    end
end
