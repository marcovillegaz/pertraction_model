classdef tUnifacGolden < matlab.unittest.TestCase
    % tUnifacGolden
    % Regression test: unifacFV.m must reproduce the pre-refactor golden
    % output exactly. If this fails, the term-extraction refactor
    % (source/calculations/thermodynamic-models/unifac/*.m) changed the
    % physics -- diff against docs/theory/thermodynamic-models/unifac.md
    % before assuming the golden value itself needs updating.
    %
    % Fixture captured from the live UNIFAC_test.m path before the
    % refactor: benzene/methylAcetate/polystyrene, T=300K, x=[0.1 0.6 0.3],
    % n=1000. See tests/fixtures/unifac_fv_golden.mat.

    properties
        CompLib
        UnifacLib
    end

    methods (TestClassSetup)
        function setupLibraries(testCase)
            testCase.CompLib = CompoundsLibrary({'benzene','methylAcetate','polystyrene'});
            testCase.UnifacLib = UNIFACLibrary();
        end
    end

    methods (Test)
        function testGoldenValueUnchanged(testCase)
            golden = load(projectPath('tests','fixtures','unifac_fv_golden.mat'));
            LnGamma = unifacFV(testCase.CompLib, testCase.UnifacLib, 300, [0.1 0.6 0.3], 1000);
            testCase.verifyEqual(LnGamma, golden.LnGamma, 'AbsTol', 1e-10);
        end

        function testFactoryMatchesDirectCall(testCase)
            % ThermoModel.create("unifac-fv", ...) must produce the same
            % result as calling unifacFV.m directly -- the class layer is
            % orchestration only, it must not change the math.
            model = ThermoModel.create("unifac-fv", testCase.CompLib, testCase.UnifacLib);
            viaModel = model.computeActivityCoefficient(300, [0.1 0.6 0.3]);
            viaFunction = unifacFV(testCase.CompLib, testCase.UnifacLib, 300, [0.1 0.6 0.3], 1000);
            testCase.verifyEqual(viaModel, viaFunction);
        end
    end
end
