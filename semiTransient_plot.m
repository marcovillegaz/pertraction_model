function [c_result,D_fick,flux] = semiTransient_plot(x_fit,S,unifac,mn_cell,Cb_0,tf,h, ...
    T,w_poly,cteEq,vol_phases,At,l,phase)
% This function employ the finite difference in the FIc's second law to
% calculate the cocnetration in both phases of a perstraction process.
% Input:
%   x_fit: row array where the the unknown variables is grouped.
%       ex: x_fit = [aij,xi,FVP,cteEq]
%   S: struct where de propeties of every compound are organized in fields
%   unifac: struct where unifac data is saved in groups 'Aij' and 'groups'
%   mn_cell: cell that contains postions of the aij unknowns in unifac.Aij
%   Cb_0: array with the molar concentration of both phases [aq ext]
%   tf: final time of modelation [h]
%   h: step time [h] (must be fraction of 1 if t_exp are integer numbers)
%   T: system temperature [K]
%   w_poly: column array with mass fraction of component in the membrane
%   cteEq: row array with the equilibria constnat of phases [eqAq eqExt]
%   vol_phases: row array with the volume of each phase [cm3] [Aq Ext]
%   At,l: membrane area an thickess [cm]
%   phase: string with the name of the phase that it concentration are 
%       going to be returned. 'aq': aqueous phase, 'ext': extractant phase
% Output:
%   c_result = column array the

%% CODE
% Put x_fit data in a object that the main function can use
[S,unifac,cteEq] = enterData(x_fit,S,unifac,mn_cell,cteEq);

%% Data extraction from structs
component = fields(S);      % compounds names
n = length(component);      % # of compounds

rho = zeros(n,1); % density [g/cm3] 
MW = zeros(n,1);  % molar wieght [g/mol]      
for i = 1:n
    rho(i) = S.(component{i}).density(T)./1000;       
    MW(i) = S.(component{i}).molarWeigth{1,4};       
end

%% Previous calculations
x_poly = (w_poly'./MW)/sum(w_poly'./MW);  % molar fraction (in the membrane)
Ct = sum(rho./MW);   % mixture molar density [mol/cm3] (in the membrane)

% Fick's diffusion coefficient using Maxwell-Stefan model
D_fick = ms_fick(S,unifac,x_poly,T);   % [cm2/s] (in the membrane)

% In the future add the diffusion coefficient in liquid phase and then
% consider the mass trasnfer resistance in all phases.

%% TIME LOOP SEMI-STATIONARY STATE
time = (0:h:tf)'.*3600;    % time vector [s]

% column(1) correspond to aq. phase, while column(2) for extractant phase
% Preallocating results
c_bulk = zeros(length(time),2);   % bulk concentration field [mol/cm3] 
flux = zeros(length(time)-1,1);     % molar flux [mol/cm2 s]  
delta_x = zeros(length(time)-1,1);   % molar gradient [-] 

% concentration at t = 0
c_bulk(1,:) = Cb_0;

for i = 1:length(time)-1 % time loop
    % Set old molar concentration
    c_bulk_old = c_bulk(i,:);
    
    % Concentration in interphase of both phases (equilibria)
    c_inter = c_bulk_old.*cteEq;

    % Transform molar concentration to mass fraction (in interphase)
    w_inter = (c_inter.*MW')./rho(end);
    w_inter(2,:) = 1 - w_inter;    % second row is polymer mass fraction

    % molar fraction in interphase (polymer side) (molar2weight)
    x_inter(:,1) = (w_inter(:,1)./MW)./sum(w_inter(:,1)./MW);
    x_inter(:,2) = (w_inter(:,2)./MW)./sum(w_inter(:,2)./MW);
    
    % Molar fracton gradient 
    delta_x(i) = x_inter(1:end-1,2) - x_inter(1:end-1,1);

    % Molar flux calculation (first Fick's law)
    flux(i) = -Ct*D_fick*(delta_x(i)./l);    % molar flux [mol/cm2 s]

    % Computing new concentration at t = t + dt
    dt = time(i+1) - time(i);

    % New concentration at aqueous (1)  and extractant (2)  phase
    c_bulk(i+1,1) = (c_bulk_old(1).*vol_phases(1) - At.*flux(i).*dt)/vol_phases(1);
    c_bulk(i+1,2) = (c_bulk_old(2).*vol_phases(2) + At.*flux(i).*dt)/vol_phases(2);
    
end

%% PLOT
% Extract concentration at the experimental times.
if strcmp(phase,"aq")
    p = 1;
elseif strcmp(phase,"ext")
    p = 2;
end

c_result = c_bulk(:,p);



%% FUNCTIONS
function D_fick = ms_fick(S,unifac,x_poly,T)
% This function employs Maxwell-Stefan diffusion therory tu calculate the
% Fick's diffusion matrix of the Fick's law. Also use the free volume
% theory to estimate the self-*difussion coefficientes. For more details,
% read the reference in the description of each function. 
% Inputs:
%   S: struct with system component data
%   unifac: struct with unifaac data 
%   x_poly: molar fraction of component inside the polymer
%   T: temepratrue of the system [K]
% Output:
%   D: Fick diffusion matrix [cm2/s]
%
% Reference of Maxwel-Stefan theory:
%   Taylor, R., & Krishna, R. (1993). Multicomponent Mass Transfer (1st ed.).
%       John Wiley & Sons.

%% CODE
% Difussion M-S matrix [B]
D0_infinite =SiddiqiLucas(T,S); % Diffusion at infinite dulution 
%fprintf('\n%-30s\n',"Diffusion coefficient at infinite dilution ... OK! ")

D_mutual = KooijmanTaylor(x_poly,D0_infinite); % Mtual difussion coefficient
%fprintf('\n%-30s\n',"Mutual diffusion coefficient for nonpolymers ... OK! ")

D_self = VrentasVrentas(S,x_poly,T); % Self diffusion 
%fprintf('\n%-30s\n',"Self diffusion coefficients ... OK! ")

D_mutual = Kubaczka(x_poly,D_self,D_mutual); % Mutual diffusion coefficient matrix (completed)
%fprintf('\n%-30s\n',"Mutual diffusion coefficient for polymers pair ... OK! ")

B = Bmatrix(D_mutual,x_poly); % B matrix calcultation
%fprintf('\n%-30s\n\n',"Diffusion coefficient matrix [B] ... OK! ")
%disp(B)

% GAMA MATRIX
Gamma = thermodynamicsFactors(S,unifac,x_poly',T);
%fprintf('%-30s\n\n',"Thermodinamic factor matrix [Gamma] ... OK! ")
%disp(Gamma)

% FICK'S DIFUSSION MATRIX [d]
D_fick = B\Gamma;   % fick difussion matrix (cm2/s)
%fprintf('%-30s\n\n',"Fick's diffusion matrix [D] ... OK! ")
%disp(D_fick)

end

function [S,unifac,cteEq] = enterData(x,S,unifac,mn_cell,cteEq)
% This function split the x cell and save the corresponding variables in
% fields S and unifac to be used in M-S diffusion model.
% Input: 
%   x: array that contain fitted variables from ga.m
%   S: struct that contains all the component and their properties
%   unifac: struct that contain the UNIFAC data needed in two fields:
%   Aij and groups
%   mn_cell: cell that contains postions of the aij unknowns in unifac.Aij
% Output:
%   S,unifac: same struct but data from x were added
%   cteEq: array with the two equilibria between bulk concentration and
%   interphase concentration in both phases. [eqAq, eqExt]

%% CODE
sz = size(mn_cell);  % # of aij unknowns
compName = fields(S);
ncomp = length(compName) -1; % # of component with unknown xi

% Splitting data
A = x(1:sz(1));
xi = x(sz(1)+1:sz(1)+ncomp);
FVP_poly = x(sz(1)+ncomp+1:sz(1)+ncomp+3);
eq = x(end);

% Putting unknown variables in unifac.Aij struct
for j = 1:sz(1)
    m = mn_cell{j,2}(1);
    k = mn_cell{j,2}(2);
    unifac.Aij(m,k) = A(j);
end

% Putting unknown variables in S.compound.FVP struct
% xi parameter
for j = 1:length(compName)-1
    S.(compName{j}).FVP{6,4} = xi(j);
end

% FV parameter fo polymer
for j = 3:5
    S.(compName{end}).FVP{j,4} = FVP_poly(j-2);
end

% Depending the case, you can unblock the lines.
% % equilibria constant in aqueous phase
% cteEq(1) = eq;  

% equilibria constant in extractant phase
cteEq(2) = eq;  % 

end

end