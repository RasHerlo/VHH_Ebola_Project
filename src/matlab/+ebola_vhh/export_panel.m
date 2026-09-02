function out = export_panel(fig, fig_id, panel_id, subpart)
%EXPORT_PANEL  Write PDF and PNG for one subpart under output/<fig>/panels/<panel>/.

    P = ebola_vhh.load_paths();
    outdir = fullfile(P.output, fig_id, 'panels', panel_id);
    if exist(outdir, 'dir') ~= 7
        mkdir(outdir);
    end

    base = fullfile(outdir, subpart);
    set(fig, 'Color', 'w');

    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
        exportgraphics(fig, [base '.png'], 'Resolution', 300);
    else
        print(fig, [base '.pdf'], '-dpdf', '-painters');
        print(fig, [base '.png'], '-dpng', '-r300');
    end

    out = struct('pdf', [base '.pdf'], 'png', [base '.png']);
end
