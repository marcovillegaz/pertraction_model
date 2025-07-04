function viscosityFunction = sastriRao(sastriRaoGroups, Tb)
% sastriRao Returns a viscosity function handle based on Sastri & Rao group
% contribution.
%
% viscosityFunction = sastriRao(groupData, Tb) creates a function handle
%   to estimate liquid viscosity [Pa·s] vs temperature T [K], using the
%   Sastri & Rao group contributions.
%
%   Inputs:
%       groupData : Nx4 matrix or table with columns:
%                  [DeltaMu, DeltaN, ocurrence]
%       Tb        : normal boiling point [K] (in otherProps Sheet)
%
%   Output:
%       viscosityFunction : function handle, so you can call
%                        mu = viscosityFunction(T);
%
%   REFERENCE:
%       Sastri, R. S., & Rao, K. K. (1992). A new group contribution method for 
%       predicting viscosity of organic liquids. 50, 9–25. 
%       https://doi.org/https://doi.org/10.1016/0300-9467(92)80002-R


%% Code
    sastriRaoGroups(:,1) = []; % Delete first column
    sastriRaoGroups = table2array(sastriRaoGroups); % Convert to numeric   
    
    % Check size
    if size(sastriRaoGroups, 2) < 3
        error('sastriRaoGroups must must have at least 3 columns [DeltaMu, DeltaN, ocurrence].');
    end
    
    % Extract group contribution coefficients
    mu_b = sum(sastriRaoGroups(:,1).*sastriRaoGroups(:,3));
    N  = 0.2 + sum(sastriRaoGroups(:,1).*sastriRaoGroups(:,3));

    % Return the function handle (verify units)
    viscosityFunction = @(T) 0.001.*mu_b.*(exp((4.5396 + 1.0309.*log(Tb)).*(1-((3-2.*(T./Tb)).^0.19./(T./Tb))-0.38.*log((T./Tb)).*(3-2.*(T./Tb)).^(-0.81)))).^(-N);

end