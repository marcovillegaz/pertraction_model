classdef UNIFACLibrary < ThermoLibrary
    properties (Access = public)
        Aij
        groups
    end
    
    methods
        %% Constructor: load unifac information into internal storage
        function obj = UNIFACLibrary(folderPath, filename)
            % Defaults inputs
            arguments
                folderPath (1,:) char = fullfile('data', 'unifac-data')
                filename (1,:) char = 'unifac.xlsx'
            end
            
            % Call parent constructor to set filePath
            obj@ThermoLibrary(folderPath, filename);
      
            % Load groups as cell array from "groups" sheet
            groupsTable = readtable(obj.filePath, 'Sheet', 'groups',"VariableNamingRule","preserve");
            obj.groups = table2cell(groupsTable);
            % Number of groups to define range for interaction parameters
            numGroups = size(obj.groups, 1);
            interactionRange = obj.getRange(numGroups);
            % Load interaction parameters matrix from "interactionParams" sheet
            obj.Aij = cell2mat(readcell(obj.filePath, 'Sheet', 'interactionParams', 'Range', interactionRange));
            
            fprintf('%-30s\n', 'UNIFAC data loaded successfully!');
        end
    end

    methods (Access = protected)
        %% getRange Generates Excel-style range for n x n matrix starting at B2
        function range = getRange(~, n)
        % getRange Generates an Excel-style cell range string for an n x n matrix.
        %   range = getRange(n) returns a string representing the Excel cell range 
        %   starting at cell B2 and covering an n-by-n matrix. For example, 
        %   getRange(3) returns 'B2:D4'.
           
            startRow = 2;  % row 2
            startCol = 2;  % column B
                        
            % Create Excel-style range dynamically
            endCol = startCol + n - 1;
            endRow = startRow + n - 1;
            % Convert numeric column to Excel letter (e.g., 2 -> 'B')
            colLetter = @(col) char('A' + col - 1);
            rangeStart = [colLetter(startCol), num2str(startRow)];
            rangeEnd   = [colLetter(endCol),   num2str(endRow)];
            range = [rangeStart ':' rangeEnd];
        end
    end
end