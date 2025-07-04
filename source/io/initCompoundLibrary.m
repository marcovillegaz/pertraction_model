function compoundLibrary = initCompoundLibrary(compound_list,folderPath)
% INITCOMPOUNDLIBRARY - Loads all compounds in \data\compound-data stored 
% as excel spreadsheet into a master struct.
%
% Input:
%   compound_list (cell) - List of compound names that defines the system
%                           (must match filenames without .xlsx)
%   folderPath (char)    - (Optional) Path to folder containing Excel files.
%                           Default is 'data/compound-data'
%
% Output:
%   compoundLibrary (struct) - Struct with fields named after compounds
%                     (e.g., compoundLibrary.Water, compoundLibrary.PCB77)

%% Default folder path
if nargin < 2 || isempty(folderPath)
    folderPath = fullfile('data', 'compound-data');
end

%% Code
compoundLibrary = struct(); % Struct definition

for k = 1:length(compound_list) % Files Loop 
    compoundName = compound_list{k};
    try
        % Full path to the Excel file
        compoundPath = fullfile(folderPath, compoundName + ".xlsx");
        % Load compound data as struct
        compound = loadCompoundData(compoundPath);
        % Generate valid fieldname
        fieldName = matlab.lang.makeValidName(compoundName); 
        % Store in master struct
        compoundLibrary.(fieldName) = compound;
        
        fprintf("✔ Loaded compound: %s\n", compoundName);
    catch ME
        warning("⚠️ Could not load %s: %s", compoundName, ME.message);
    end
end

fprintf("✅ Loaded %d compounds into compoundLibrary.\n", length(fieldnames(compoundLibrary)));
end
