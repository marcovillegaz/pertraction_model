function [fitFunc, stats] = propertyFit(dataTable, fitmodel, T_range,varargin)
% Fits experimental data (property vs. temperature) from a table using nonlinear regression.
%
% Input:
%   dataTbl  : Table with at least three columns: T_C, T_K, property
%   fitmodel : A string specifying the fitting model, e.g., 'exp2', 'poly2', etc.
%   T_range  : 2-element array with temperature range in Kelvin, e.g., [273 300]
%
% Optional Inputs (name-value pairs):
%   'PropertyName' : A string used for plot labels. Default: 'Property'.
%   'Plot'         : Boolean, whether to display plot. Default: true.
% 
% Output:
%   fitFunc  : Function handle for the fitted model
%   coeffs   : Coefficients of the fitted model
%   r2       : Coefficient of determination of the fit
%
% ADVICES: 
% THis function does not perform unit conversion

%% Default
    %Default values for optional inputs 
    propertyName = 'Property';
    doPlot = false;

    % Parse optional inputs
    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'propertyname'
                propertyName = varargin{i+1};
            case 'plot'
                doPlot = varargin{i+1};
            case 'compoundname'
                compoundName = varargin{i+1};
            otherwise
                error('Unknown parameter: %s', varargin{i});
        end
    end

%% Code
    
    % Keep only the first three columns (T_C, T_K, Property)
    dataTable = dataTable(:, 1:3);
    % Rename for consistency
    dataTable.Properties.VariableNames = {'T_C', 'T_K', 'Property'};
    % Remove rows with any missing or non-numeric values
    dataTable = rmmissing(dataTable);
    % Filter temperature range
    isInRange = dataTable.T_K >= T_range(1) & dataTable.T_K <= T_range(2);
    dataTable = dataTable(isInRange, :);
    
    % Check if there's enough data
    if height(dataTable) < 3
        error("Not enough data points in the selected temperature range.");
    end

    % Extract from table
    T_K = dataTable.T_K;
    y = dataTable.Property;

    % Perform fitting
    [fitObj, stats] = fit(T_K, y, fitmodel);
    
    % Cfit into function handle
    fitFunc = @(x) feval(fitObj, x);
    
    % save fit results
    saveFitResults(fitObj, stats, compoundName, propertyName)

    % plot and save (optional)
    if doPlot
        plotAndSaveFit(T_K, y, fitFunc, compoundName, propertyName)
    end
end

%% Private functions
function saveFitResults(fitObj, stats, compoundName, propertyName)
%SAVEFITRESULTS Save fit object and stats to a .txt file

    % Ensure output folder exists
    outputFolder = fullfile('data','fit-results');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Clean filename parts
    compoundName = lower(compoundName);
    propertyName = lower(propertyName);

    % Generate file name
    txtFileName = fullfile(outputFolder, ...
        strcat(compoundName, '_', propertyName, '_fit.txt'));

    % Capture text outputs
    fitSummary = evalc('disp(fitObj)');
    statsSummary = evalc('disp(stats)');

    % Write to file
    fid = fopen(txtFileName, 'w');
    if fid == -1
        warning('Could not open file for writing: %s', txtFileName);
        return;
    end

    fprintf(fid, 'Fitted Model:\n%s\n', fitSummary);
    fprintf(fid, 'Fit Statistics:\n%s\n', statsSummary);
    fclose(fid);

    fprintf('\nFit results were save in :%s', txtFileName);
end

function plotAndSaveFit(T_K, y, fitFunc, compoundName, propertyName)
%PLOTANDSAVEFIT Create and save a plot of the fitted data (invisible figure)

    % Ensure the output folder exists
    outputFolder = 'images';
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Lowercase and clean names
    compoundName = lower(compoundName);
    propertyName = lower(propertyName);

    % Create invisible figure
    f = figure('Visible', 'off');
    ax = axes(f);

    % Plot data
    plot(ax, T_K, y, 'o', 'DisplayName', 'Data');
    hold(ax, 'on');

    % Smooth curve vector
    T_fit = linspace(min(T_K), max(T_K), 200);
    plot(ax, T_fit, fitFunc(T_fit), '-', 'DisplayName', 'Fit');
    hold(ax, 'off');

    % Labels
    xlabel(ax, 'Temperature [K]');
    switch propertyName
        case 'density'
            ylabel(ax, 'Density [kg/m^3]');
        case 'viscosity'
            ylabel(ax, 'Viscosity [Pa·s]');
        otherwise
            ylabel(ax, propertyName);
    end

    legend(ax, 'Location', 'best');

    % File name
    fileName = fullfile(outputFolder, strcat(compoundName, '_', propertyName, '_fit.jpg'));

    % Export and close
    exportgraphics(f, fileName, 'ContentType', 'image', 'Resolution', 300);
    close(f);

    fprintf('\n%s vs temperature plot was save in: %s', propertyName,fileName)
end