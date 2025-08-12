classdef ThermoLibrary
    % ThermoLibrary Base class for thermodynamic data libraries
    %   Handles common logic for loading data files from a folder.
    
    properties (Access = public)
        filePath 
    end
    
    methods
        function obj = ThermoLibrary(folderPath, filename)
            % Constructor: set folder path, using default if not provided
            arguments
                folderPath (1,:) char
                filename (1,:) char = ''
            end
            obj.filePath = fullfile(folderPath, filename);
            fprintf('%-30s %s\n', 'ThermoLibrary loading data from:', obj.filePath);
        end
    end
    
end