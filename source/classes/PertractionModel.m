% This class employs the finite difference in the Fick's second law to
% calculate the concetration in both phases of a perstraction process
% equation (X). 

% Class definition
classdef PertractionModel
     properties
        Temperature           % [K]
        MembraneThickness     % [m]
        MembraneArea          % [m^2]
        Volumes               % struct with fields: V_aq, V_org, etc. [m^3]
        Compounds             % Array of Compound objects
        UNIFACModel           % UNIFACModel object
        PolymerComposition    % Mass fractions or mole fractions of polymer groups
    end

    methods
        function obj = PerstractionModel(T, L, A, volumes, compounds, unifacModel, wPoly)
            % Constructor
            obj.Temperature = T;
            obj.MembraneThickness = L;
            obj.MembraneArea = A;
            obj.Volumes = volumes;
            obj.Compounds = compounds;
            obj.UNIFACModel = unifacModel;
            obj.PolymerComposition = wPoly;
        end

        function D = getDiffusionCoefficients(obj)
            % Placeholder method – could call a diffusion model
            % (e.g., Vrentas-Duda, Siddiqi-Lucas)
            D = zeros(length(obj.Compounds), 1);
            for i = 1:length(obj.Compounds)
                D(i) = obj.Compounds(i).estimateDiffusivity(obj.Temperature, obj.PolymerComposition);
            end
        end

        function gamma = getActivityCoefficients(obj, x)
            % Calls UNIFAC model
            gamma = obj.UNIFACModel.computeActivityCoeffs(x, obj.Temperature);
        end

        function obj = updateFromFitVector(obj, x_fit)
            % This could parse x_fit to update diffusivities, Aij, etc.
            % e.g., x_fit = [D1 D2 ... Aij1 Aij2 ...]
            % Then assign to respective compound/UNIFAC parameters
            obj.UNIFACModel = obj.UNIFACModel.updateAij(x_fit);
            for i = 1:length(obj.Compounds)
                obj.Compounds(i) = obj.Compounds(i).updateFromFitVector(x_fit);
            end
        end
    end
end

