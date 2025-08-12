classdef ThermoModel
    % ThermoModel
    % Base class for thermodynamic models.
    % Stores references to thermodynamic and compound libraries.
    properties (Access = public)
        thermoLib     % Library with thermodynamic parameters (e.g., UNIFAC groups)
        compoundsLib  % Library with compound information (e.g., molecular data)
    end 

    methods
        function obj = ThermoModel(compoundsLibObj, thermoLibObj)
            arguments
                compoundsLibObj {mustBeA(compoundsLibObj, 'CompoundsLibrary')}
                thermoLibObj   {mustBeA(thermoLibObj, 'ThermoLibrary')}
            end

            obj.compoundsLib = compoundsLibObj;            
            obj.thermoLib = thermoLibObj;
        end
        
        function Lngamma = computeActivityCoefficient(~, varargin)
            % computeActivityCoefficient
            % General interface for computing activity coefficients.
            % Child classes will override this method with model-specific implementations.
            
            error('ThermoModel:NotImplemented', ...
                'computeActivityCoefficient must be implemented in the child class.');
        end
    end
end 