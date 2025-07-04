function UNIFAC = loadUnifacData(filename, folderPath, dataRange)
%LOADUNIFACDATA Loads UNIFAC group interaction data from Excel.
%
%   UNIFAC = loadUnifacData(filename)
%   UNIFAC = loadUnifacData(filename, folderPath)
%   UNIFAC = loadUnifacData(filename, folderPath, dataRange)
%
%   Inputs:
%     - filename:   Name of the Excel file to read (e.g., 'groups.xlsx')
%     - folderPath: Path to the folder containing the file (default: 'data/unifac-data')
%     - dataRange:  Excel range to extract group data (default: 'A3:O8')
%
%   Output:
%     - UNIFAC: Struct with fields:
%         - Aij: Matrix of group interaction parameters
%         - groups: Cell array of group metadata (IDs, names, etc.)

%% Default
% Set default range if not provided
if nargin < 2 || isempty(dataRange)
    dataRange = 'A3:O8';
end

% Set default folder path if not provided
if nargin < 3 || isempty(folderPath)
    folderPath = fullfile('data', 'unifac-data');
end


%% Load UNIFAC Data
filePath = fullfile(folderPath, filename);
fprintf('%-30s %s\n', 'Extracting UNIFAC group data from:', filePath);

% ENHENCE THIS PROCESS IN THE FUTURE TO BE MORE GENERAL
% Create empty struct
UNIFAC = struct();  
% Extract group data
UNIFAC.groups = table2cell(readtable(filePath,'Sheet', "groups"))

% get number of groups
[numGroups,~] = size(UNIFAC.groups)
% generate the range in interactionParams sheet
intractionParamsRange = getRange(numGroups)

% extract interaction parameters matrix
UNIFAC.Aij = cell2mat(readcell(filePath, 'Range', intractionParamsRange, 'Sheet', "interactionParams"))

fprintf('%-30s\n', "OK!");

end

function range = getRange(n)
% getRange Generates an Excel-style cell range string for an n x n matrix.
%   range = getRange(n) returns a string representing the Excel cell range 
%   starting at cell B2 and covering an n-by-n matrix. For example, 
%   getRange(3) returns 'B2:D4'.

    % Convert n to a cell range string (e.g., 'B2:F6')
    startRow = 2;  % starting from row 2
    startCol = 2;  % 'B' column = column 2
    % Create Excel-style range dynamically
    endCol = startCol + n - 1;
    endRow = startRow + n - 1;
    % Convert numeric column to Excel letter (e.g., 2 -> 'B')
    colLetter = @(col) char('A' + col - 1);
    rangeStart = [colLetter(startCol), num2str(startRow)];
    rangeEnd   = [colLetter(endCol),   num2str(endRow)];
    range = [rangeStart ':' rangeEnd]
end 