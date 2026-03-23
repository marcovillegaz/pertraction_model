classdef Membrane
    properties
        area double              % m²
        thickness double         % m
        temperature double       % K
        molarFraction double     % Composition inside membrane
        diffusivityModel         % Object of type DiffusivityModel
    end

    methods
        function obj = initialize(obj, w_poly)
            obj.molarFraction = computeXpoly(w_poly);
            obj.Ct = computeMolarDensity();
            obj.D_fick = obj.computeDiffusivity();
        end

        function D = computeDiffusivity(obj)
            % encapsulate ms_fick() here
        end
    end
end
