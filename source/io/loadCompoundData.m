function compound = loadCompoundData(filePath)
% LOADCOMPOUNDDATA - Loads compound data from a structured Excel file.
%
% Syntax:
%   compound = loadCompoundData(compoundName)
%
% Description:
%   Loads the sheets 'OtherProps', 'Density', and 'Viscosity' (if available)
%   from a compound-specific Excel file located in 'data/compound-data/'.
%   If Density or Viscosity are missing, it will estimate them using
%   critical properties from OtherProps. If critical properties are also missing,
%   an error message is shown.
%
% Inputs:
%   filePath - Full path to the Excel file (e.g., "test_compounds/Water.xlsx")
%
% Outputs:
%   compound - Struct containing fields:
%     .name           - Compound name
%     .OtherProps     - Table of critical or general properties
%     .Density        - Table (if exists) or estimated table (if estimated)
%     .Viscosity      - Table (if exists) or estimated table (if estimated)
%
% Example:
%   water = loadCompoundData("Water");


    %% CONSTANTS
    ranges = {'A1:F20','A3:O8'};   % The second range is for UNIFAC
    
    %%
    if ~isfile(filePath)
        error("File not found: %s", filePath);
    end
    
    % Extract compound name from the file name
    [~, compoundName, ~] = fileparts(filePath);

    % Get list of sheet in compounds.xlsx
    [~, sheetNames] = xlsfinfo(filePath);
    
    % Define .name field
    compound.name = compoundName; 

    %% Load OtherProps (required for estimation)
    if ismember("OtherProps", sheetNames)
        % Extract aand rearrange data
        % compound.OtherProps = readtable(filePath, "Sheet", "OtherProps");
        otherPropsCell = readcell(filePath,'Range',ranges{1},'sheet',"OtherProps")
        [otherPropsCell,varGroups] = rearrangeCell(otherPropsCell);

        % Save otherProps in compound struct
        for i = 2:length(varGroups)
            GroupArray = otherPropsCell{i};
            % GroupArray(:,pos2cut) = [] ;
            compound.(varGroups{i}) = GroupArray ; % save in struct
        end

    else
        warning("OtherProps sheet missing for compound %s", compoundName);
        compound.OtherProps = table();  % Empty table
    end

    %% Load or estimate Density
    if ismember("Density", sheetNames)
        % Read experimental density
        densityTable = readtable(filePath, "Sheet", "Density");

        % Fit experimental data 
        [fitFunc, ~] = propertyFit(densityTable, "exp2",[283.15 363.15],...
            'compoundName',compound.name,'propertyName',"Density","Plot",true);
 
        % Save density ad funciton handle
        compound.density = fitFunc;
 
    else
        % If there is not experimental data, density is etsimated
        fprintf("No experimental Density data for %s. Attempting to estimate...\n", compoundName);
        
        if ismember("GCVOL60",sheetNames) 
            elbroGroups = readtable(filePath,'Sheet',"GCVOL60");
            compound.density = GCVOL60(elbroGroups, compound.molarWeigth);

        else
            error("Cannot estimate Density for %s: missing GCVOL60 sheet", compoundName);
        end
    end

    %% Load or estimate Viscosity
    if ismember("Viscosity", sheetNames)
        % Read experimental density
        viscosityTable = readtable(filePath, "Sheet", "Viscosity");
        % Fit experimental data 
        [fitFunc, ~] = propertyFit(viscosityTable, "exp2",[283.15 363.15],...
            'compoundName',compound.name,'propertyName',"Viscosity","Plot",true);

        % Save density ad funciton handle
        compound.Viscosity = fitFunc;

    else
        fprintf("No experimental Viscosity data for %s. Attempting to estimate...\n", compoundName);
        if ismember("SastriRao",sheetNames)
            sastriRaoGroups = readtable(filePath,'Sheet',"SastriRao");
            Tb = compound.boilingPoint{1,4};
            compound.Viscosity = sastriRao(sastriRaoGroups, Tb);
        else
            error("Cannot estimate Viscosity for %s: missing critical properties.", compoundName);
        end
    end


end

%% Private functions
function [groupedData, groupNames] = rearrangeCell(propCell)
%REARRANGECELL Reorganizes property cell array into grouped sub-cells.
%
%   [groupedData, groupNames] = rearrangeCell(propCell) takes a cell array
%   from the 'OtherProps' sheet and organizes it by property groups, making
%   it easier to access and process each group individually.
%
%   Inputs:
%       propCell    - Cell array where each row represents a property and the
%                     first column indicates group names (e.g., 'molarWeight').
%                     Subsequent rows contain property data, with <missing> in
%                     the first column where the group is continued.
%
%   Outputs:
%       groupedData - Cell array of sub-cells, each containing properties of a
%                     specific group, excluding the group name column.
%       groupNames  - Cell array of group names (non-missing values from first column).

    fprintf('%-30s', "Reorganizing property data...")

    % Size of the input cell array
    [numRows, ~] = size(propCell);

    % Extract the first column containing group labels
    groupColumn = propCell(:,1);

    % Find rows where group names are NOT missing
    isGroupRow = ~cellfun(@(x) isa(x, 'missing'), groupColumn);
    groupIndices = find(isGroupRow);

    % Extract group names
    groupNames = groupColumn(groupIndices);

    % Add one more index to allow slicing up to the end of the last group
    groupIndices(end+1) = numRows;

    % Preallocate the output cell
    groupedData = cell(length(groupNames), 1);

    % Loop through group indices and extract corresponding blocks (excluding group column)
    for i = 1:length(groupNames)
        startRow = groupIndices(i);
        if i<length(groupNames)
            endRow   = groupIndices(i+1) - 2; % exclude the next group row and trailing <missing>
            groupedData{i} = propCell(startRow:endRow, 2:end); % exclude the first column
        else  %Last group
            groupedData{i} = propCell(startRow:end, 2:end); % exclude the first column
        end
    end

    fprintf('%-10s\n', "Done!")
end






