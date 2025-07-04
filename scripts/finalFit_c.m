% This script employs genetic algorithm to fit experimental perstraction
% concentrations to a mass transfer model that uses Maxwel-Stefan diffusion
% theory and Vrentas-Vrentas free volume theory. After the fit, the results
% are saved. 

clc, clear
%% SYSTEM INPUT VARIABLES
T = 40+273.15;      % system temperature [K]
tf = 7;             % final time of experiment [h]    
h = 0.25;           % time step (must be fraction of 1 if t_exp are integer numbers)                      
At = 27.2484;       % mass transfer area [cm2]
l = 0.009;          % membrane thickness [cm]
MW_pcb = 291.99;     % [g/mol]  density of PCB77

vol_phases = [75 75];          % volume of phases [cm3]
rho_phases = [0.9563 1.320];   % density of phases [g/cm3] 

% cell with positin of unknowns parameters in Aij matrix
mn_cell = {      
    'a_43_25',[6 5]
    'a_25_42',[5 6]
    };

%% LOADING NECESSARY DATA AND PRE-PROCESSING
% Load equilibrium data
load eqData.mat       % eqData w_poly
    cteEq = [eqData.POMS 0];     % Equilibria constant of aq phase. [mol/cm3] / [mol/cm3]

% Load data stored in structs
load compoundDataFVP.mat   % S unifac
    unifac = UNIFAC_data;   
    S.poms.FVP{6,4} = 1;        %% TAKE A LOOK
    S.poms.density = @(T) 910;  %% TAKE A LOOK

% Loas experimental data stored in struct
load expPerstract.mat  % expPerstract
    t_exp = expPerstract.time;   % experimental times
    w.aq = (expPerstract.waterPOMS.aq).*10^-6;   % mass fractions
    w.ext = (expPerstract.waterPOMS.ext).*10^-6;
    c_exp.aq = (w.aq.*rho_phases(1))./MW_pcb;    % molar concentrations
    c_exp.ext = (w.ext.*rho_phases(2))./MW_pcb;
    Cb_0 = [c_exp.aq(1,1) c_exp.ext(1,1)];  % Molar concentration at t = 0 for both phases [aq ext]

%% LINEAR CONSTRAINS AND BOUNDS FOR GENETIC ALGORITHM
% unknowns group interaction parameters of UNIFAC
aij_lb = zeros(1,2) - 10000
aij_ub = zeros(1,2) + 10000

% xi parameter of the free volume theory of Vrentas
xi_lb = [0];   
xi_ub = [1];

% Free volume parameters of the free volume theory of Vrentas
FVP_poly_lb = [0 0 -500]  
FVP_poly_ub = [1e-10 1e-3 500]

% Unknown equilibria constant of phase
eq_lb = [0]
eq_ub = [1000] 

% algortihm constrains
A = [];
b = [];
Aeq = [];
beq = [];

lb = [aij_lb,xi_lb,FVP_poly_lb,eq_lb]
ub = [aij_ub,xi_ub,FVP_poly_ub,eq_ub]

nonlcon = [];

constrains = {A,b,Aeq,beq,lb,ub,nonlcon};

%% GENETIC ALGORITHM OPTIONS (optional)
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
    'MaxTime',300, ...                % input
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

% Generating option varType
options = optimoptions("ga",displayOptions{:}, ...
stopOptions{:},funcOptions{:},inputOptions{:});


    %% Model to fit experimental data
    % Second fick's law is modelated using Maxwel-Stefan diffusion for a
    % perstraction experiment in semi-transient state. 
    % semiTransient.m emmployes a time loop to preict the molar concentration 
    % in both side of the membrane.
    
    ymodel = @semiTransient_fit_both;   
    ydata = [c_exp.aq(:,1);c_exp.ext(:,1)]  % important
    xdata = t_exp
    
    % Minimal residual squered: function to optimize
    fun = @(cte) sum((ydata-ymodel(cte,S,unifac,mn_cell,Cb_0,t_exp,tf,h, ...
        T,w_poly,cteEq,vol_phases,At,l)).^2); 

    %% Genetic algorithm
    [cte,~] = ga(fun,7,constrains{:},options);
    close

    %% Coefficient of determination
    r2 = 1-sum((ydata-ymodel(cte,S,unifac,mn_cell,Cb_0,t_exp,tf,h, ...
        T,w_poly,cteEq,vol_phases,At,l)).^2)/sum((ydata-mean(ydata)).^2);

%% RESULTS PROCESSING AND SAVE
cteNames = {'a_43_25';'a_25_42';'xip';'D0_i';'K_1i/gamma';
    'K_2i - T_gi';'eqExt';'r2'};
cteUnits = {'K';'K';'-';'cm2/s';'cm3/g K';'K';'-';'-'};
cteCell = num2cell([cte';r2]);
colName = {'variable','unit','value'};
cteFinal = [cteNames,cteUnits,cteCell];
cteFinal = [colName;cteFinal];

% Save results in .txt file
writecell(cteFinal,'waterPOMS_fit_c.txt')
clear
    