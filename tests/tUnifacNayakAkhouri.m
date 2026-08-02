classdef tUnifacNayakAkhouri < matlab.unittest.TestCase
    % tUnifacNayakAkhouri
    % Placeholder literature-validation test against
    % data/input/unifac/test_nayakAkhouri.xlsx (group/interaction-parameter
    % table, apparently a Nayak & Akhouri UNIFAC-FV validation case --
    % unreferenced elsewhere in the codebase). Currently only checks that
    % the table loads with the expected shape.
    %
    % TODO: once the literature reference and its reported activity
    % coefficients are available, replace testTableLoads with a real
    % regression test: build a CompoundsLibrary matching the paper's
    % system, call the appropriate unifac*.m variant, and verifyEqual
    % against the paper's reported values (with a documented tolerance).

    methods (Test)
        function testTableLoads(testCase)
            unifacLib = UNIFACLibrary(projectPath('data','input','unifac'), 'test_nayakAkhouri.xlsx');
            testCase.verifyGreaterThan(size(unifacLib.groups,1), 0, ...
                'expected at least one group row');
            testCase.verifyEqual(size(unifacLib.Aij,1), size(unifacLib.Aij,2), ...
                'interaction parameter matrix should be square');
        end
    end
end
