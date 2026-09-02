function setup()
%SETUP  Add this repo's MATLAB package to the path and ensure local config exists.
%
%   Run once per MATLAB session from the repository root:
%       setup

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src', 'matlab'));

    cfg = fullfile(root, 'config', 'paths.json');
    example = fullfile(root, 'config', 'paths.example.json');
    if exist(cfg, 'file') ~= 2
        copyfile(example, cfg);
        fprintf('Created config/paths.json from the example.\n');
    end

    out = fullfile(root, 'output');
    if exist(out, 'dir') ~= 7
        mkdir(out);
    end

    fprintf('VHH Ebola figures ready. Try:\n');
    fprintf('  ebola_vhh.list_panels\n');
    fprintf('  ebola_vhh.fig5.E.cfd\n');
    fprintf('  ebola_vhh.fig5.assemble\n');
end
