function params = unifacGroupParams(compoundsLibrary, unifacData, T, n)
%UNIFACGROUPPARAMS Shared UNIFAC group-parameter setup.
%   params = unifacGroupParams(compoundsLibrary, unifacData, T, n) extracts
%   and rearranges the group-contribution data (R, Q, group occurrences,
%   interaction parameters) and derives the per-compound volume/area
%   parameters (r, q, l) shared by the combinatorial, residual, and
%   free-volume terms (unifacCombinatorial.m, unifacResidual.m,
%   unifacFreeVolume.m).
%
%   Inputs:
%     compoundsLibrary - CompoundsLibrary object
%     unifacData        - UNIFACLibrary object (fields .Aij, .groups)
%     T                  - absolute temperature [K]
%     n                  - degree of polymerization applied to the LAST
%                          compound (the polymer), following
%                          r(end) = n*r(end). Pass n = 1 (the default) for
%                          systems with no polymer component -- this
%                          leaves r unchanged, matching legacy UNIFAC.
%
%   Output:
%     params - struct with fields g, c, z, r, q, l, Q, V, psi, v
%
% Reference: Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).
%   The Properties of Gases and Liquids (5th ed.). McGraw-Hill, ch. 8.10.

arguments
    compoundsLibrary
    unifacData
    T (1,1) double
    n (1,1) double = 1
end

%% Extract data
Aij = unifacData.Aij;            % group interaction parameters
groupData = unifacData.groups;   % cell with groups data

MW = compoundsLibrary.extractPropertyAsArray({"molarWeigth",1});
rho_fh = compoundsLibrary.extractPropertyAsFunction("density");

groupData(:,6) = [];             % deleting MW column
[g,c] = size(groupData);         % g = # of groups
c = c-5;                         % c = # of components

v = MW./rho_fh(T);               % molar volume vector [cm3/mol]

z = 10;                           % UNIFAC coordination number

%% R, Q, and l
r = zeros(c,1);
q = zeros(c,1);
l = zeros(c,1);
aux = cell2mat(groupData(:,4:end));

for k = 1:c
    r(k) = sum(aux(:,1).*aux(:,k+2));
    q(k) = sum(aux(:,2).*aux(:,k+2));
    l(k) = (z/2)*(r(k)-q(k))-(r(k)-1);
end

% Polymer volume correction (no-op when n = 1)
r(end) = n*r(end);

%% Residual-term group data (rows = groups, columns = compounds)
Q = cell2mat(groupData(:,5));
V = cell2mat(groupData(:,6:end));
psi = exp(-Aij./T);

params = struct('g',g, 'c',c, 'z',z, 'r',r, 'q',q, 'l',l, 'Q',Q, 'V',V, 'psi',psi, 'v',v);
end
