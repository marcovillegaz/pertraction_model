function LnGamma = activityCoefficients(modelType, compoundsLibrary, unifacData, x, T, n)
%ACTIVITYCOEFFICIENT Compute activity coefficients of the compounds in mixture 
% using a specified model and return a array. 
%   LnGamma = activityCoefficient(S, x, T, n, 'Model', 'vdw-fv', 'Data', UNIFAC_data)

% Parse name-value inputs
    switch lower(modelType)
        case 'unifac-vdw-fv'
            LnGamma = UNIFAC_vdW_FV(compoundsLibrary, unifacData, x, T, n);
        case 'unifac-test'
            LnGamma = UNIFAC_test(compoundsLibrary, unifacData, x, T, n);
        otherwise
            error('Unsupported UNIFAC model type: %s', modelType);
    end
end