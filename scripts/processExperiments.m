% This script process the experimental data stored in
% experimental_final.xlsx to be used in the matlab codes. 

% The data is then used as input in the mass transfer model and to cretae
% plots with experimental points. 

clc, clear
fileName = 'final_experiments.xlsx';
filePath = fullfile(pwd, 'data', 'input','experimental', fileName);

% Name of the sheet wich contians the relevant experimental information. 
EXPERIMENT_LIST = {'milk_POMS','water_PEBA','milk_PEBA','water_POMS'};

% Name of the sheet that contains the experimental data to obtain the
% partition constant.
PARTITION_SHEET = "Partition Constants";

%% Initialize
fprintf('PROCESESSING %12s\n\n',fileName)

%% PROCESS DATA FROM PERTRACTION EXPERIMENTS
% Initialize struct
experimentalStruct = struct();

% Loop over each experimental run (sheet)
for s = 1:length(EXPERIMENT_LIST)
    sheet = EXPERIMENT_LIST{s};
    fprintf('[%s] (%d of %d) Processing experiment: %s\n', ...
        datestr(now, 'HH:MM:SS'), s, length(EXPERIMENT_LIST), sheet);
    rawData = readtable(filePath, 'Sheet', sheet);
    
    % Clean and prepare
    rawData = rmmissing(rawData);  % Remove missing
    rawData.PHASE = string(rawData.PHASE);
    
    % Unique phases
    phases = unique(rawData.PHASE);
    
    % Loop by phase
    for i = 1:length(phases)
        phase = phases(i);
        phaseData = rawData(rawData.PHASE == phase, :);

        % Group by time
        [groupedTimes, ~, groupIdx] = unique(phaseData.TIME);
        meanConc = accumarray(groupIdx, phaseData.CONCENTRATION, [], @mean);
        stdConc = accumarray(groupIdx, phaseData.CONCENTRATION, [], @std);

        % Sort by time
        [sortedTime, sortIdx] = sort(groupedTimes);
        meanConc = meanConc(sortIdx);
        stdConc = stdConc(sortIdx);

        % Save in struct
        experimentalStruct.(sheet).(phase).time = sortedTime;
        experimentalStruct.(sheet).(phase).mean_concentration = meanConc;
        experimentalStruct.(sheet).(phase).std_concentration = stdConc;
    end
end

save(fullfile(pwd,'data','experimentalStruct.mat'),'experimentalStruct')

%% PROCESS DATA FROM PARTITION CONSTANT EXPERIMENTS
fprintf('[%s] Processing experiment: %s\n', datestr(now, 'HH:MM:SS'),PARTITION_SHEET);

table = readtable(filePath, 'Sheet', PARTITION_SHEET);
table.Properties.VariableNames = matlab.lang.makeValidName(table.Properties.VariableNames);

% Compute group mean using groupsummary (like Pandas groupby().mean())
summaryTable = groupsummary(table, {'MEMBRANE', 'PHASE', 'TIME'}, {'mean', 'std'}, 'CONCENTRATION_ug_mL_');
% Rename columns for clarity
summaryTable.Properties.VariableNames{'mean_CONCENTRATION_ug_mL_'} = 'MEAN_CONCENTRATION';
summaryTable.Properties.VariableNames{'std_CONCENTRATION_ug_mL_'} = 'STD_CONCENTRATION';

% Load necessary constants
Vliq = 3;       % liquid volume in mL
Vmem = 0.009;   % membrane volume in mL 
MW = 291.99;     % molar weight g/mol of PCB77

% run("system_constants.m")

partitionTable = computePartition(summaryTable, Vliq, Vmem,PCB77_MW);
disp(partitionTable)

save(fullfile(pwd,'data','partitionConstants.mat'),'partitionTable')

fprintf('PROCESESSING COMPLETE! %12s\n\n',fileName)
clear

%% LOCAL FUNCTIONS
function resultTable = computePartition(meanTable, Vliq, Vmem,MW)

    % Initialize result storage
    result = [];
    % Get list of membranes
    membranes = unique(meanTable.MEMBRANE);
    % Get list of each liquid phase (assuming water and extractant are both liquid)
    phases = unique(meanTable.PHASE);

    for i = 1:numel(membranes) % Membrane loop
        mem = membranes{i};

        for j = 1:numel(phases) % Phase loop
            phase = phases{j};

            % Filter data for this membrane and phase
            subset = meanTable(strcmp(meanTable.MEMBRANE, mem) & ...
                               strcmp(meanTable.PHASE, phase), :);

            % Get concentration at time 0 and 24
            Concentration0hrs = subset.MEAN_CONCENTRATION(subset.TIME == 0);
            Concentration24hrs = subset.MEAN_CONCENTRATION(subset.TIME == 24);
            
            % In case of empty cell
            if isempty(Concentration0hrs) || isempty(Concentration24hrs)
                warning('Missing time points for %s - %s. Skipping.', mem, phase);
                continue;
            end
            
            % Revisar undiades esta en kmol/m3 en el paper

            % Computations
            absorbesMass = (Concentration0hrs - Concentration24hrs)*Vliq;   % Mass absorbed (g)
            ConcentrationInMembrane = ((absorbesMass*1e-6)/MW) / Vmem;      % Concentration in the membrane (mol/mL)
            Concentration24hrs = (Concentration24hrs*1e-6)/MW ;             % g/mL to mol/ml
            % Compute partition constant
            Keq = ConcentrationInMembrane / Concentration24hrs;

            % Store result
            result = [result; {mem, phase, Keq}];
        end
    end

    % Convert to table
    resultTable = cell2table(result, ...
        'VariableNames', {'MEMBRANE', 'PHASE', 'Keq'});
end


