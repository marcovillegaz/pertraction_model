% This script plot the estimated and experimental values of the three fit 
% aq, both and ext.

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

tf = 7

% LOAD FITTED DATA
bothData = readcell('waterPOMS_fit_c.txt');
Fit.both = cell2mat(bothData(2:end-1,3))

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

% Load equilibrium data
load eqData.mat       % eqData w_poly
    cteEq = [eqData.POMS 0];     % Equilibria constant of aq phase. [mol/cm3] / [mol/cm3]

%% Estimatin cocnentration by model
t_model = [0:0.1:tf]'
fields = fieldnames(Fit)
names = {'fit_c.jpg'}

for i = 1:length(names)
    figure 
    hold on,grid on 
    xlabel('Tiempo h','FontSize',10)
    ylabel('Concecentración ug/mL','FontSize',10)

    % plotting experimetnal data
    errorbar(t_exp,expPerstract.waterPOMS.ext(:,1), ...
        expPerstract.waterPOMS.ext(:,2),"s", ...
        'MarkerFaceColor',"#D95319", ...
        'MarkerEdgeColor',"#D95319", ...
        'Color',"#D95319")
    
    errorbar(t_exp,expPerstract.waterPOMS.aq(:,1), ...
       expPerstract.waterPOMS.aq(:,2),"^", ...
        'MarkerFaceColor',"#0072BD", ...
        'MarkerEdgeColor',"#0072BD", ...
        'Color',"#0072BD")
    
    % COmpute estimated data
    cte = Fit.(fields{i})
    [c_model,D_fick_aq(i),flux] = semiTransient_plot_both(cte,S,unifac,mn_cell,Cb_0,tf,0.1, ...
        T,w_poly,cteEq,vol_phases,At,l)

    % molar cocnentration to ppm
    cppm.model.aq = (c_model(:,1).*S.pcb77.molarWeigth{1,4}./rho_phases(1))*10^6
    cppm.model.ext = (c_model(:,2).*S.pcb77.molarWeigth{1,4}./rho_phases(2))*10^6

    % plotting estimated data by model
    plot(t_model,cppm.model.aq,'-', ...
        'Color',"#0072BD",'LineWidth',1)
    plot(t_model,cppm.model.ext,'-', ...
        'Color',"#D95319",'LineWidth',1)

    legend('Fase extractante','Fase acuosa', ...
        'FontSize',10)

    % save image
    exportgraphics(gca,names{i},'Resolution',800)
    close
end

