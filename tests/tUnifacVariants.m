classdef tUnifacVariants < matlab.unittest.TestCase
    % tUnifacVariants
    % Per-variant smoke and shape tests for the three UNIFAC models
    % (legacy, FV, vdW-FV), exercised through ThermoModel.create so the
    % factory dispatch itself is covered too.

    properties
        CompLib
        UnifacLib
        T = 300
        X = [0.1 0.6 0.3]
    end

    methods (TestClassSetup)
        function setupLibraries(testCase)
            testCase.CompLib = CompoundsLibrary({'benzene','methylAcetate','polystyrene'});
            testCase.UnifacLib = UNIFACLibrary();
        end
    end

    methods (Test)
        function testLegacyShapeAndFiniteness(testCase)
            g = testCase.computeVia("unifac");
            testCase.verifySize(g, [3 1]);
            testCase.verifyTrue(all(isfinite(g)), 'legacy output must be finite');
            testCase.verifyTrue(~any(isnan(g)), 'legacy output must not contain NaN');
        end

        function testFVShapeAndFiniteness(testCase)
            g = testCase.computeVia("unifac-fv");
            testCase.verifySize(g, [3 1]);
            testCase.verifyTrue(all(isfinite(g)), 'FV output must be finite');
            testCase.verifyTrue(~any(isnan(g)), 'FV output must not contain NaN');
        end

        function testVdwFVShapeAndFiniteness(testCase)
            g = testCase.computeVia("unifac-vdw-fv");
            testCase.verifySize(g, [3 1]);
            testCase.verifyTrue(all(isfinite(g)), 'vdW-FV output must be finite');
            testCase.verifyTrue(~any(isnan(g)), 'vdW-FV output must not contain NaN');
        end

        function testLegacyDiffersFromFV(testCase)
            % The free-volume term must actually change the answer --
            % otherwise unifacFreeVolume.m isn't being applied.
            gLegacy = testCase.computeVia("unifac");
            gFV = testCase.computeVia("unifac-fv");
            testCase.verifyTrue(any(abs(gLegacy - gFV) > 1e-6), ...
                'legacy and FV should not produce identical output');
        end

        function testFactoryReturnsExpectedClass(testCase)
            testCase.verifyClass(ThermoModel.create("unifac", testCase.CompLib, testCase.UnifacLib), ...
                'UNIFACLegacyModel');
            testCase.verifyClass(ThermoModel.create("unifac-fv", testCase.CompLib, testCase.UnifacLib), ...
                'UNIFACFVModel');
            testCase.verifyClass(ThermoModel.create("unifac-vdw-fv", testCase.CompLib, testCase.UnifacLib), ...
                'UNIFACvdWFVModel');
        end

        function testUnknownVariantErrors(testCase)
            testCase.verifyError(...
                @() ThermoModel.create("not-a-real-model", testCase.CompLib, testCase.UnifacLib), ...
                'ThermoModel:UnknownVariant');
        end
    end

    methods (Access = private)
        function g = computeVia(testCase, variant)
            model = ThermoModel.create(variant, testCase.CompLib, testCase.UnifacLib);
            g = model.computeActivityCoefficient(testCase.T, testCase.X);
        end
    end
end
