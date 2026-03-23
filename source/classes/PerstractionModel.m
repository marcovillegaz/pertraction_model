classdef PerstractionModel
    properties
        compoundsLibrary   % struct or object array with compound info
        unifacLibrary      % struct with UNIFAC groups and parameters
        systemConfig       % struct: T, V1, V2, A, L, polymer_mass_fraction, etc.
        equilibriumConstants   % vector, fitted or default

        % Derived
        x_poly             % molar fraction of polymer phase
        Ct_poly            % total molar concentration in polymer
    end

    methods
        function obj = PerstractionModel(compoundsLibrary, unifacLibrary, systemConfig)
            obj.compoundsLibrary = compoundsLibrary;
            obj.unifacLibrary = unifacLibrary;
            obj.systemConfig = systemConfig;
        end

        function obj = applyFit(obj, x_fit)
            % Unpack fit vector into K, UNIFAC Aij, free volume, etc.
            % You might want to refactor this if you want multiple models.
            % For now, store equilibrium constants and update libraries.
            obj.equilibriumConstants = x_fit.K;

            % Update unifacLibrary, compoundLibrary based on x_fit
            obj.unifacLibrary = updateUnifac(obj.unifacLibrary, x_fit.Aij);
            obj.compoundsLibrary = updateCompounds(obj.compoundsLibrary, x_fit.FV_params);

            % Recompute mixture properties
            obj = obj.computeMixtureProperties();
        end

        function obj = computeMixtureProperties(obj)
            % Calculates x_poly and Ct_poly from compound info and polymer mass fraction
            T = obj.systemConfig.T;
            w_poly = obj.systemConfig.polymer_mass_fraction;

            MW = arrayfun(@(c) c.MW, obj.compoundsLibrary);
            rho = arrayfun(@(c) c.getDensity(T), obj.compoundsLibrary);

            mass_total = (1 - w_poly) ./ rho + w_poly / obj.systemConfig.polymer_density;
            mol_total  = (1 - w_poly) ./ MW;
            obj.Ct_poly = mol_total ./ mass_total;

            obj.x_poly = mol_total / sum(mol_total);
        end

        function K = getPartitionCoefficients(obj)
            % Could be calculated or just returned from stored values
            K = obj.equilibriumConstants;
        end

        function [c_poly_eq] = getEquilibriumConcentrations(obj, c_bulk)
            % Applies partition coefficients to bulk concentrations
            K = obj.getPartitionCoefficients();
            c_poly_eq = K .* c_bulk;  % element-wise
        end
    end
end
