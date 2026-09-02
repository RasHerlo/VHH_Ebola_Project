function P = load_paths()
%LOAD_PATHS  Read config/paths.json and add repo-relative output paths.

    repo = ebola_vhh.repo_root();
    cfg = fullfile(repo, 'config', 'paths.json');
    if exist(cfg, 'file') ~= 2
        example = fullfile(repo, 'config', 'paths.example.json');
        copyfile(example, cfg);
    end

    P = jsondecode(fileread(cfg));
    P.repo_root = repo;
    P.output = fullfile(repo, 'output');
    P.manifest = fullfile(repo, 'config', 'figure_manifest.json');
end
