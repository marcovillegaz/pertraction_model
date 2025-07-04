% Author: Marco Villegas
% Contact: marco.villegas@usach.cl

clc,clear
%% Input parameters
T = (0:5:60)'+ 273.15;          % [K] Temperature range
nga = 5;  % number of times that the genetic algortihm is applied
plotNames = {'PCB77','Water','[omim][Tf2N]','Acetonitrile'};
datafile = 'compoundData.mat';

%% Linear constrains and bounds
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0 -500];   % lower limits of FVP based in bibliography
ub = [1e-2 1e-2 500];  % upper limits of FVP based in bibliography
nonlcon = [];

constrains = {A,b,Aeq,beq,lb,ub,nonlcon};

%% GENETIC ALGORTIHM OPTIONS (optional)
inputOptions = {
    'PopulationSize',200, ...        % input
    'EliteCount',10, ... 
    'CrossoverFraction',0.8
    };
funcOptions = {
    'CreationFcn',@gacreationuniform, ...  
    'SelectionFcn',@selectionroulette, ... 
    'FitnessScalingFcn',@fitscalingrank, ... 
    'CrossoverFcn',@crossoverlaplace, ... 
    'MutationFcn',@mutationadaptfeasible, ... 
    'HybridFcn',@fmincon
    };

stopOptions = {
    'MaxGenerations',Inf, ...        % input
    'MaxTime',0.5, ...                % input
    'FitnessLimit',1e-20, ... 
    'MaxStallGenerations',50, ... 
    'MaxStallTime',Inf, ... 
    'FunctionTolerance',1e-100, ... 
    'ConstraintTolerance',1e-3
    };

displayOptions = {
     'PlotFcn',[], ... 
     'Display','iter'
     };

%% MAIN CODE
fprintf('%30s\n',repmat('=',7))
fprintf('%30s\n','FV PARAMETERS FITTING')
fprintf('%30s\n',repmat('-',7))

%% loading data
fprintf('%-30s',strcat("Importing compound data from ",string(datafile)))
load(datafile);
component = fields(S);
component(end) = [];
fprintf('%-30s\n'," OK!")
fprintf('%30s\n',repmat('-',7))

%% Component loop
best_results = cell(length(component),4); % preallocation 

for i = 1:length(component)
    % Extract data froms atruct
    fprintf('%-30s',strcat("Extracting data of ",string(component{i})," ..."))
    
    MW = S.(component{i}).molarWeigth{1,4};  % [g/mol]   molar weight
    Vc = S.(component{i}).criticalProp{3,4}; % [cm3/mol] molar volume of solvent at critical temperature
    V0 = S.(component{i}).FVP{2,4};          % [cm3/g]   specific volume of liquid solvent at 0 K - Haward, R. N. (1970). Table 3
    
    fprintf('%-30s\n'," OK!")

    % Generating data to fit
    fprintf('%-30s',"Generating data to fit ... ")

    if i == 1 
        % only for the case of PCB77, where density function is not 
        % vectorized because the nature fo GCVOL60 model
        rho = zeros(length(T),1);
        mu = zeros(length(T),1);
        for j = 1:length(T)
            rho(j) = S.(component{i}).density(T(j))./1000;  % [g/cm3] density array
            mu(j) = S.(component{i}).viscosity(T(j));      % [Pa s] viscosity array
        end

    else
        rho = S.(component{i}).density(T)./1000;  % [g/cm3] density array
        mu = S.(component{i}).viscosity(T);      % [Pa s] viscosity array
    end

    fprintf('%-30s\n'," OK!")

    %% Free volume parameters model
    const = 0.124e-16;  % [mol^(2/3)]  Vrentas model constant
    R = 8.314e6;        % [cm3 Pa/mol K]  Ideal gas constant*
    n = 3;              % Number of model constant (Free volume parameters)
    
    % Model to fit
    ymodel = @(cte,xdata) cte(1).*(exp(-V0./(cte(2).*(cte(3)+xdata))));
    
    % Variable change
    xdata = T;
    ydata = (const.*(Vc.^(2/3)).*rho.*R.*xdata)./(MW.*mu);
    
    % Residual sum of squares
    fun = @(cte) sum((ydata-ymodel(cte,xdata)).^2); % TARGET FUNCTION

    % Generating option varType for Genetic Algorithm
    options = optimoptions("ga");
    options = optimoptions(options,displayOptions{:}, ...
        stopOptions{:},funcOptions{:},inputOptions{:});
    
    %% GA loop
    fprintf('\n%-30s\n',strcat("STARTING GA LOOP FOR: ",string(component{i})))

    results = cell(nga,4);
    for j = 1:nga
        % Genetic algorithm
        [cte,~] = ga(fun,3,constrains{:},options);close
        % Coefficient of determination
        r2 = 1-sum((ydata-ymodel(cte,xdata)).^2)/sum((ydata-mean(ydata)).^2);
        
        fprintf('%-30s\n',strcat("iteration ",num2str(j),"   r2:",num2str(r2)))
        r2 = num2cell(r2);
        cte = num2cell(cte);

        % save results
        results(j,1:3) = cte;
        results(j,4) = r2;
    end
    
    % Find best fit
    [~,idx] = max(cell2mat(results(:,4)));  
    best_results(i,:) = results(idx,:);

    %% Save FV parameters in struct
    fprintf('\n%30s',"Saving best result in struct ...")

    S.(component{i}).FVP{3,4} = best_results{i,1};
    S.(component{i}).FVP{4,4} = best_results{i,2};
    S.(component{i}).FVP{5,4} = best_results{i,3};

    fprintf('%-30s\n'," OK!")

    %% Result plot
    hold on
    xName = '$\frac{{0.124 \times {{10}^{ - 16}}V_{Ci}^{2/3}RT{\rho _i}\left( T \right)}}{{{M_1}{\eta _i}\left( T \right)}}$';
    xlabel('Temperature [K]')
    ylabel(xName,'interpreter','latex','fontsize',18)

    times = linspace(xdata(1),xdata(end));
    cte = cell2mat(best_results(i,1:3));
    plot(xdata,ydata,'ko')
    plot(times,ymodel(cte,times),'b-')

    legend('Experimental points','GA optimization','Location','northwest')
    
    imagename = strcat(component{i},"_FVP",".jpg");

    fprintf('\n%-30s',strcat("saving plot as ",imagename," ..."))
        exportgraphics(gca,imagename,'Resolution',800)

    fprintf('%-30s\n\n'," OK!")
    fprintf('%30s\n',repmat('-',7))

    close
end

save compoundDataFVP.mat S propGroups UNIFAC_data

fprintf('\n%30s',"Saving FVP bests results in .txt file")
    colnames = {'component','D0','K_1i/gamma','K_2i - T_gi','r2'};
    best_results = [component,best_results];
    best_results = [colnames;best_results];
    writecell(best_results,'FVP_results.txt')
    
fprintf('%-30s\n'," OK!")
fprintf('%30s\n',repmat('-',7))
fprintf('%30s\n','DONE')
fprintf('%30s\n',repmat('=',7))

