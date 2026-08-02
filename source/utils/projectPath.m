function p = projectPath(varargin)
%PROJECTPATH Resolve a path relative to the repository root.
%   p = projectPath('data','input','compounds') returns an absolute path
%   to data/input/compounds, independent of the current working directory.
%
%   The repo root is resolved from this file's own location
%   (source/utils/projectPath.m -> two levels up), so callers no longer
%   need to be run with CWD == repo root.
%
%   Usage:
%       compoundsFolder = projectPath('data','input','compounds');
%       unifacFile = projectPath('data','input','unifac','unifac-test.xlsx');

    thisFile = mfilename('fullpath');          % .../source/utils/projectPath
    utilsDir = fileparts(thisFile);            % .../source/utils
    sourceDir = fileparts(utilsDir);           % .../source
    root = fileparts(sourceDir);               % repo root

    p = fullfile(root, varargin{:});
end
