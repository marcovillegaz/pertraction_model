function [S,varGroups,UNIFAC] = propertiesImporter(filename,sheets,ranges)
% This funciton import all the propeties saved in the a excel file that
% follow a certain format. For more details open the file 'Properties.xlsx'
% attached to this code.
% Input
%   fileName: name of the excel spreadsheet
%   sheets: name of the sheets of the spreadhsheet
%   ranges: ranges where the data is stored in each sheet
% Output:
%   S: struct that contains the compounds as fields and the group
%   of propeties. Each compounds has the propertiees gruops as fields.
%   componentNames: cell with the component names
%   varGroupss: cell with the name of the different group of variables
%   UNIFAC:  struct that contains the UNIFAC data of the system
%
% Example of component struct
%   S is a struct with fields 
%
%    ETHANOL: [1×1 struct]
%      WATER: [1×1 struct]
%       PDMS: [1×1 struct]
%
%   where S.WATER is a struct with fields:
% 
%         molarWeigth: {2×5 cell}
%        criticalProp: {3×5 cell}
%     glassTransition: {1×5 cell}
%                 FVP: {7×5 cell}
%
%   each of this fields are a group of variables.
%
% 
% Example of UNIFAC struct
%   UNIFAC is a struct with fields 
%
%       Aij: {4×6 cell}           group interaction parameters matrix
%    groups: {4×9 cell}           groups parameters R,Q,MW and v_i (see UNIFAC.m)   

%% CODE
fprintf('%30s\n',repmat('=',7))
fprintf('%30s\n','PROPERTY MATRIX IMPORTER')
fprintf('%30s\n',repmat('-',7))

fprintf('%-30s\n\n',strcat("importing system data from ",string(filename)))

%% VARIABLE GROUPS DATA
fprintf('%-30s\n',strcat("extracting variables grouped from ",string(sheets{1})," sheet"))
propData = readcell(filename,'Range',ranges{1},'sheet',sheets{1})

%% Rearrangment of propData
fprintf('%-30s',"rearrangement...")
sz = size(propData)

% name of components
componentNames = lower(propData(1,5:end-1)') ;
n_component = length(componentNames) ; % # of components

% name of variables groups
aux = propData(:,1);
mask = cellfun(@(x) any(isa(x,'missing')), propData(:,1));
idx = find(mask == 0);
varGroups = aux(idx);

% cut groups
idx(end+1) = sz(1)+2;
c = cell(length(varGroups),1);
for i = 1:length(idx)-1
    c{i} = propData(idx(i):(idx(i+1)-2),2:end);
end

% Save and cut data
flag = 1:n_component;
flag = flag + 3;  % column where numeric data is saved
for i = 1:length(componentNames)
    for j = 1:length(varGroups)
        aux = c{j};
        pos2cut = flag  ;
        pos2cut(i) = [] ; % column to be deleted
        aux(:,pos2cut) = [] ;
        S.(componentNames{i}).(varGroups{j}) = aux ; % save in struct
    end
end
fprintf('%-30s\n\n'," OK!")



%% UNIFAC DATA
fprintf('%-30s\n',strcat("extracting UNIFAC group data from ",string(sheets{2})," sheet"))
unifacData = readcell(filename,'Range',ranges{2},'sheet',sheets{2}); 

%% Rearrangment of unifacData
fprintf('%-30s',"rearrangement...")

[n_groups,m] = size(unifacData);         % n = # of groups

UNIFAC.Aij = cell2mat(unifacData(:,(m-n_groups+1):end));  % interaction parameter between groups
UNIFAC.groups = unifacData(:,1:(m-n_groups-1));  % cell that contains UNIFAC group data
fprintf('%-30s\n'," OK!")


fprintf('%30s\n',repmat('-',7))
fprintf('%30s\n','DONE')
fprintf('%30s\n',repmat('=',7))
end


