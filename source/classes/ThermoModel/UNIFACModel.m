classdef UNIFACModel < ThermoModel
    % UNIFACModel
    % Implements activity coefficient calculation using the UNIFAC method.
    
    methods
        %% Constructor
        function obj = UNIFACModel(compoundsLibObj, thermoLibObj)
            % Call the parent constructor
            obj@ThermoModel(compoundsLibObj, thermoLibObj);
        end
        %% Compute natural logarithm of activity coefficient
        function LnGamma = computeActivityCoefficient(obj, ...
                temperature, molarWeight)
            
            % COMPUTE USING UNIFAC CALCULATIONS
            LnGamma = UNIFAC_test(obj.compoundsLib,obj.thermoLib,temperature,molarWeight,1000);
        end
    end
end
