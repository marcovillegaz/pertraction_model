classdef CompoundsLibrary
    %COMPOUNDSLIBRARY Load and manage compound data from Excel files
    %
    %   This class loads compound property data from a folder of Excel files 
    %   into a structured library. Each file should contain the data for a 
    %   single compound. Compounds are accessed by name, and data is stored 
    %   in a private struct.
    %
    %   Usage:
    %       lib = CompoundsLibrary(compoundList)
    %       lib = CompoundsLibrary(compoundList, folderPath)
    %
    %   Inputs:
    %       compoundList - Cell array of compound names (file names without .xlsx)
    %       folderPath   - Optional path to folder containing Excel files
    %
    %   Public Methods:
    %       list()       - Return names of all loaded compounds
    %       get(name)    - Return struct of a compound by name
    %
    %   Example:
    %       list = {'benzene', 'toluene'};
    %       lib = CompoundsLibrary(list, 'data/my-compounds');
    %       benzene = lib.get('benzene');

    properties (Access = private)
        compounds  % Struct with all compound data
    end

    methods
        %% Constructor: load compounds into internal storage
        function obj = CompoundsLibrary(compoundList, folderPath)
            % Defaults inputs
            arguments
                compoundList (1,:) cell
                folderPath (1,:) char = projectPath('data', 'input', 'compounds')
            end

            % Struct definition
            data = struct();
            % Files loop 
            for k = 1:length(compoundList)
                name = compoundList{k};
                try
                    % Full path to the Excel file
                    path = fullfile(folderPath, name + ".xlsx");
                    % Load compound data as struct
                    compound = loadCompoundData(path);    % ./source/io
                    % Generate valid fieldname
                    fieldName = matlab.lang.makeValidName(name);
                    % Store in data struct
                    data.(fieldName) = compound;

                    fprintf("✔ Loaded compound: %s\n", name);
                catch ME
                    warning("⚠️ Could not load %s: %s", name, ME.message);
                end
            end
            
            % define compounds struct
            obj.compounds = data;
            fprintf("✅ Loaded %d compounds into CompoundsLibrary.\n", length(fieldnames(data)));
        end

        %% Return all compound names in library
        function names = list(obj)        
            names = fieldnames(obj.compounds);
        end

        %% Return compound struct by name (case-insensitive)
        function compound = get(obj, name)            
            fieldName = matlab.lang.makeValidName(name);
            if isfield(obj.compounds, fieldName)
                compound = obj.compounds.(fieldName);
            else
                error("Compound '%s' not found in library.", name);
            end
        end
        
        %% Extract properties array based in property path 
        function propertyArray = extractPropertyAsArray(obj, propertyPath)
            % Extracts a specified property from a list of compounsd in the 
            % compoundLibrary and return a cell or array depending on the data type. 
            %
            % Inputs:
            %   obj. 
            %   propertyPath    - cell array describing the path to the property,
            %           For cell fields, propertyPath = {'fieldName', row}
            %           For function handles, propertyPath = {'fieldName'}
            % Output:
            %   propertyArray - cell array or numeric array of the extracted values
            
            compoundNames = obj.list();
            n = length(compoundNames);
            values = cell(n,1);  % Use cell in case the property is non-numeric or varies
            
            for i = 1:n
                try
                    compoundData = obj.get(compoundNames{i});
                    fieldData = compoundData.(propertyPath{1});
                    
                    % Differentiate between cell and function handle
                    if iscell(fieldData) && length(propertyPath) == 2
                        values{i} = fieldData{propertyPath{2}, 4};
                    else
                        values{i} = fieldData;
                    end
            
                catch
                    values{i} = NaN;  % Assign NaN if any issue arises
                end
            end
            
            % Try to convert to numeric if all entries are numeric
            if all(cellfun(@isnumeric, values))
                propertyArray = cell2mat(values);
            end

        end
        
        function propertyFunc = extractPropertyAsFunction(obj,propertyName)
            %EXTRACTPROPERTYASFUNCTION Returns a function handle that evaluates a property across compounds.
            %   fh = extractPropertyAsFunction(compoundLibrary, propertyName)
            %
            %   Inputs:
            %     - obj.
            %     - propertyName: Name of the property (e.g., 'density') which is a function handle in each compound.
            %
            %   Output:
            %     - fh: A function handle that takes T as input and returns a vector:
            %           [compound1.property(T), compound2.property(T), ...]
            
                compoundNames = obj.list();  % e.g., {'compound1', 'compound2'}
                n = numel(compoundNames);
            
                % Preload all function handles for efficiency
                fh_vector = cell(n, 1);
                for i = 1:n
                    compound = obj.get(compoundNames{i});
                    if isfield(compound, propertyName)
                        fh_vector{i} = compound.(propertyName);
                    else
                        fh_vector{i} = @(T) NaN;  % Fallback if field is missing
                    end
                end
                
                disp(fh_vector)
                
                % Return a unified function that evaluates all
                propertyFunc = @(T) cellfun(@(f) f(T), fh_vector);
            end 
    end
end

