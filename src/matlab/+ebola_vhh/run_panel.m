function run_panel(fig_id, panel_id)
%RUN_PANEL  Execute every subpart of one panel, then compose it.
%
%   ebola_vhh.run_panel('fig5', 'E')
%   ebola_vhh.run_panel('s3', 'C')

    M = ebola_vhh.load_manifest();
    if ~isfield(M.figures, fig_id)
        error('Unknown figure id: %s', fig_id);
    end
    F = M.figures.(fig_id);
    if ~isfield(F.panels, panel_id)
        error('Unknown panel %s in %s', panel_id, fig_id);
    end
    Pan = F.panels.(panel_id);
    parts = ebola_vhh.as_cellstr(Pan.subparts);

    for i = 1:numel(parts)
        fn = str2func(['ebola_vhh.' fig_id '.' panel_id '.' parts{i}]);
        fprintf('Running %s.%s.%s ...\n', fig_id, panel_id, parts{i});
        fn();
    end

    fprintf('Composing %s %s ...\n', fig_id, panel_id);
    ebola_vhh.compose_panel(fig_id, panel_id);
end
