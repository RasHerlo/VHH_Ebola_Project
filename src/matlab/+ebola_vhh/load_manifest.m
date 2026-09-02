function M = load_manifest()
%LOAD_MANIFEST  Figure/panel map from config/figure_manifest.json.

    P = ebola_vhh.load_paths();
    M = jsondecode(fileread(P.manifest));
end
