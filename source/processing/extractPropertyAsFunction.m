function fh = extractPropertyAsFunction(compoundLibrary, propertyName)
%EXTRACTPROPERTYASFUNCTION Returns a function handle that evaluates a property across compounds.
%
%   fh = extractPropertyAsFunction(compoundLibrary, propertyName)
%
%   Inputs:
%     - compoundLibrary: Struct with compound names as fields.
%     - propertyName: Name of the property (e.g., 'density') which is a function handle in each compound.
%
%   Output:
%     - fh: A function handle that takes T as input and returns a vector:
%           [compound1.property(T), compound2.property(T), ...]

    compoundNames = fieldnames(compoundLibrary);  % e.g., {'compound1', 'compound2'}
    n = numel(compoundNames);

    % Preload all function handles for efficiency
    fh_vector = cell(n, 1);
    for i = 1:n
        compound = compoundLibrary.(compoundNames{i});
        if isfield(compound, propertyName)
            fh_vector{i} = compound.(propertyName);
        else
            fh_vector{i} = @(T) NaN;  % Fallback if field is missing
        end
    end
    
    disp(fh_vector)
    
    % Return a unified function that evaluates all
    fh = @(T) cellfun(@(f) f(T), fh_vector);
end