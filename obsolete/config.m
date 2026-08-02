function cfg = config()
    % Singleton pattern for configuration
    persistent cachedCfg
    
    if isempty(cachedCfg)
        % Only executed the first time
        cachedCfg = struct();

        % Define your settings here
        cachedCfg.polymerList = {'polystyrene', 'polyethylene', 'PVC'};
        cachedCfg.dataFolder = fullfile('data', 'compound-data');
        cachedCfg.resultsFolder = fullfile('results');

        % You could even load something from file if needed
        % cachedCfg.someData = readtable('config_values.xlsx');
    end

    % Return the cached instance
    cfg = cachedCfg;
end