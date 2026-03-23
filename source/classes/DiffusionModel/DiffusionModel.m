classdef DiffusionModel
    %DIFFUSIONMODEL Computes diffusion coefficients using a defined method
    % This class interfaces with a ThermodynamicLibrary or CompoundLibrary
    % to calculate Maxwell-Stefan or other models.

    properties
        compLib % object handle to CompoundLibrary
    end

    methods
        function obj = DiffusionModel(compoundLibrary)
            obj.compLib = compoundLibrary;
        end
    end
    
    methods (Abstract)
        D = computeDiffusivity(obj, compoundA, compoundB, T, x)
    end

end