function S = TynCalus(S)
% This function estimates the liquid molar volume at the normal boiling
% point using the Tyn & Calus(1975) correlation.
% Input:
%   component: struct that contains 
% Output:
%   componen: cell with where the volume at normal boiling point (cm3/mol)
%   has been added

%% Formula 
fprintf('%30s\n',repmat('=',7))
fprintf('%30s\n','TYN & CALUS MOLAR VOLUME')
fprintf('%30s\n',repmat('-',7))

fprintf('%-30s\n',"calculating molar volume at Tb ... OK!")

component = fieldnames(S);   % field names
for i = 1:length(component)-1
    Vc = S.(component{i}).criticalProp{3,4};
    Vb = 0.285*Vc^1.048;
    S.(component{i}).boilingPoint{2,4} = Vb;
end

fprintf('%30s\n',repmat('-',7))
fprintf('%30s\n','DONE')
fprintf('%30s\n',repmat('=',7))

end 