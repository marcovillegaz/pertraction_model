function values = extractPropertyAsArray(compoundLibrary, propertyPath)
% Extracts a specified property from a list of compounsd in the 
% compoundLibrary and return a cell or array depending on the data type. 
%
% Inputs:
%   compoundLibrary - struct of compounds
%   compoundNames   - cell array of compound names (strings)
%   propertyPath    - cell array describing the path to the property,
%           For cell fields, propertyPath = {'fieldName', row}
%           For function handles, propertyPath = {'fieldName'}
%
% Output:
%   values - cell array or numeric array of the extracted values

%% CODE
compoundNames = fieldnames(compoundLibrary);
n = length(compoundNames);
values = cell(n,1);  % Use cell in case the property is non-numeric or varies

for i = 1:n
    try
        compoundData = compoundLibrary.(compoundNames{i});
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
    values = cell2mat(values);
end

end