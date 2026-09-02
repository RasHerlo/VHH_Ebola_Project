function root = repo_root()
%REPO_ROOT  Absolute path to the git repository root.

    here = fileparts(mfilename('fullpath'));
    % .../src/matlab/+ebola_vhh
    root = fileparts(fileparts(fileparts(here)));
end
