function results = runAllTests()
%RUNALLTESTS Run the full test suite for the pertraction model.
%   results = runAllTests() adds source/ to the MATLAB path (anchored to
%   this file's own location, so it works from any current working
%   directory) and then runs every test in this folder via the built-in
%   matlab.unittest framework.
%
%   Usage:
%       cd tests; runAllTests
%   or, from anywhere:
%       run(fullfile(repoRoot,'tests','runAllTests.m'))
%
%   When called with no output argument (e.g. from `matlab -batch`), a
%   test failure raises an error so the batch process exits non-zero --
%   useful for CI. When called as `r = runAllTests()`, failures are
%   returned in r without throwing, for interactive inspection.

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
addpath(genpath(fullfile(repoRoot,'source')));
addpath(thisDir);

results = runtests(thisDir);
disp(results);

if nargout == 0 && any([results.Failed])
    error('runAllTests:TestsFailed', '%d of %d tests failed.', ...
        sum([results.Failed]), numel(results));
end
end
