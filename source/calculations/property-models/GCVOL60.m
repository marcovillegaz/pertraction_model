function densityFunction = GCVOL60(elbroGroups, MW)
% GCVOL60 Returns a density function handle based on Elbro group contribution.
%
%   densityFunction = GCVOL60(elbroGroups, MW) creates a function handle
%   to estimate the liquid density [kg/m^3] as a function of temperature T [K],
%   using group contribution data from Elbro et al. (1991).
%
%   Input:
%       elbroGroups : Table that contain elbro groups of the compound
%       MW          : Molecular weight of the compound [g/mol]
%
%   Output:
%       densityFunction : function handle that calculates density(T) [kg/m^3]
%
%   Usage:
%       f = GCVOL60(elbroGroups, 146.14);
%       rho = f(298.15);  % get density at 25 °C
%
%   Reference:
%       Ihmels, E. C. (2003). Extension and Revision of the Group Contribution 
%       Method GCVOL for the Prediction of Pure Compound Liquid Densities.
%       Ind. Eng. Chem. Res., 42(2), 408–412.
%       https://doi.org/10.1021/IE020492J
%% Code
    elbroGroups(:,1) = []; % Delete first column
    elbroGroups = table2array(elbroGroups); % Convert to numeric
    
    % Check size
    if size(elbroGroups, 2) ~= 4
        error('elbroGroups must be an Nx4 matrix with columns [A, B, C, n].');
    end

    % Extract group contribution coefficients
    A = elbroGroups(:,1);
    B = elbroGroups(:,2);
    C = elbroGroups(:,3);
    n = elbroGroups(:,4);

    % Return the function handle
    densityFunction = @(T) 1000 .* (MW ./ sum((A + B .* T + C .* T.^2) .* n));

end