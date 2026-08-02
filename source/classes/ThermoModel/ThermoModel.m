classdef ThermoModel
    % ThermoModel
    % Abstract base class for thermodynamic (activity coefficient) models.
    % Stores references to thermodynamic and compound libraries, and
    % provides ThermoModel.create(...) to construct a named model variant
    % without the caller needing to know the concrete subclass.
    %
    % Subclasses must implement:
    %   computeActivityCoefficient(obj, T, x) -> LnGamma column vector
    %   RequiredCompoundProps (Constant) -> string array of compound
    %       struct field names this model needs (empty if none)

    properties (Access = public)
        thermoLib     % Library with thermodynamic parameters (e.g., UNIFAC groups)
        compoundsLib  % Library with compound information (e.g., molecular data)
    end

    properties (Abstract, Constant)
        RequiredCompoundProps  % string array; declares what this model needs from each compound
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

        function validateCompoundData(obj)
            % VALIDATECOMPOUNDDATA Check that every compound in
            %   compoundsLib carries the fields this model declares via
            %   RequiredCompoundProps, and that required FVP rows are
            %   present and non-empty. Raises a clear error naming the
            %   compound and the missing property -- deliberately does
            %   NOT use CompoundsLibrary.extractPropertyAsArray, whose
            %   bare catch (CompoundsLibrary.m) silently turns a missing
            %   field into NaN instead of erroring.
            required = obj.RequiredCompoundProps;
            if isempty(required)
                return
            end

            names = obj.compoundsLib.list();
            for i = 1:numel(names)
                compound = obj.compoundsLib.get(names{i});
                for p = 1:numel(required)
                    propName = required(p);
                    if ~isfield(compound, propName)
                        error('ThermoModel:MissingCompoundData', ...
                            'Compound "%s" is missing required property "%s" for model %s.', ...
                            names{i}, propName, class(obj));
                    end

                    if propName == "FVP"
                        % FVP is a 7x5 cell (see
                        % docs/theory/diffusion-models/free-volume-theory.md).
                        % Rows 2 (V*), 3 (D0_i), 4 (K1/gamma), 5
                        % (K2-Tg), 7 (E_i) must have a value in column 4.
                        % Row 6 (xi_ip) is known-empty in the current
                        % data and is intentionally excluded.
                        fvp = compound.FVP;
                        requiredRows = [2 3 4 5 7];
                        for r = requiredRows
                            cellValue = fvp{r,4};
                            if isempty(cellValue) || ismissing(cellValue)
                                error('ThermoModel:MissingCompoundData', ...
                                    'Compound "%s" has an empty FVP row %d (column 4) required by model %s.', ...
                                    names{i}, r, class(obj));
                            end
                        end
                    end
                end
            end
        end
    end

    methods (Abstract)
        LnGamma = computeActivityCoefficient(obj, temperature, molarFraction)
    end

    methods (Static)
        function obj = create(variantName, compoundsLibObj, thermoLibObj, opts)
            % THERMOMODEL.CREATE Factory for UNIFAC model variants.
            %   model = ThermoModel.create("unifac-fv", compLib, unifacLib)
            %   model = ThermoModel.create("unifac-fv", compLib, unifacLib, ...
            %                               "PolymerizationDegree", 500)
            %
            %   Valid variantName values: "unifac" (legacy, no free-volume
            %   term), "unifac-fv", "unifac-vdw-fv".
            %
            %   Validates the compound library against the chosen
            %   variant's requirements before returning, so a missing
            %   FVP/density fails fast at construction rather than deep
            %   inside a later calculation.
            arguments
                variantName (1,1) string
                compoundsLibObj {mustBeA(compoundsLibObj, 'CompoundsLibrary')}
                thermoLibObj   {mustBeA(thermoLibObj, 'ThermoLibrary')}
                opts.PolymerizationDegree (1,1) double = 1000
            end

            switch lower(variantName)
                case {"unifac", "unifac-legacy", "legacy"}
                    obj = UNIFACLegacyModel(compoundsLibObj, thermoLibObj);
                case {"unifac-fv", "fv"}
                    obj = UNIFACFVModel(compoundsLibObj, thermoLibObj, opts.PolymerizationDegree);
                case {"unifac-vdw-fv", "vdw-fv"}
                    obj = UNIFACvdWFVModel(compoundsLibObj, thermoLibObj, opts.PolymerizationDegree);
                otherwise
                    error('ThermoModel:UnknownVariant', ...
                        ['Unknown thermodynamic model variant "%s". ' ...
                         'Valid options: "unifac", "unifac-fv", "unifac-vdw-fv".'], variantName);
            end

            obj.validateCompoundData();
        end
    end
end
