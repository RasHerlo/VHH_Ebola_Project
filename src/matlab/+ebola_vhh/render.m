function fig = render(path_key, fig_id, panel_id, subpart)
%RENDER  Open a configured source file and export it as a named subpart.
%
%   Baseline helper: current panels start as the 2022 .fig / TIFF files.
%   Replace the body of a panel function with new plotting code when that
%   subpart is being redrawn.

    P = ebola_vhh.load_paths();
    if ~isfield(P.files, path_key)
        error('ebola_vhh:unknownKey', 'No files.%s in config/paths.json', path_key);
    end
    fig = ebola_vhh.open_source(P.files.(path_key));
    ebola_vhh.export_panel(fig, fig_id, panel_id, subpart);
end
