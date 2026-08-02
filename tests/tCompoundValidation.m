classdef tCompoundValidation < matlab.unittest.TestCase
    % tCompoundValidation
    % Verifies ThermoModel.create fails fast, with a clear error, when a
    % compound library doesn't carry the data a variant needs -- rather
    % than silently propagating NaN into a downstream calculation.
    %
    % Uses real compound files rather than synthetic fixtures, since two
    % of them already exercise the two failure modes that matter:
    %   - omimTf2N.xlsx has NO "FVP" group at all in its OtherProps sheet.
    %   - Water.xlsx HAS an "FVP" group, but rows 3 (D0_i), 4 (K1/gamma),
    %     5 (K2-Tg) are empty -- "present but empty" is a distinct bug
    %     class from "field missing" and both must be caught.

    properties
        UnifacLib
    end

    methods (TestClassSetup)
        function setupLibraries(testCase)
            testCase.UnifacLib = UNIFACLibrary();
        end
    end

    methods (Test)
        function testMissingFVPFieldRejectedByFV(testCase)
            compLib = CompoundsLibrary({'omimTf2N'}, projectPath('data','input','compounds'));
            testCase.verifyError(...
                @() ThermoModel.create("unifac-fv", compLib, testCase.UnifacLib), ...
                'ThermoModel:MissingCompoundData');
        end

        function testMissingFVPFieldAcceptedByLegacy(testCase)
            % Legacy UNIFAC declares no FVP requirement, so the same
            % compound library that fails FV must succeed here.
            compLib = CompoundsLibrary({'omimTf2N'}, projectPath('data','input','compounds'));
            model = ThermoModel.create("unifac", compLib, testCase.UnifacLib);
            testCase.verifyClass(model, 'UNIFACLegacyModel');
        end

        function testEmptyFVPRowsRejectedByFV(testCase)
            compLib = CompoundsLibrary({'Water'}, projectPath('data','input','compounds'));
            testCase.verifyError(...
                @() ThermoModel.create("unifac-fv", compLib, testCase.UnifacLib), ...
                'ThermoModel:MissingCompoundData');
        end

        function testFullyPopulatedCompoundsAcceptedByFV(testCase)
            % Sanity check the positive case with the compounds actually
            % used by main.m -- all three carry a fully populated FVP.
            compLib = CompoundsLibrary({'benzene','methylAcetate','polystyrene'});
            model = ThermoModel.create("unifac-fv", compLib, testCase.UnifacLib);
            testCase.verifyClass(model, 'UNIFACFVModel');
        end
    end
end
