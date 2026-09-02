function fig = assemble_figure(fig_id, varargin)
%ASSEMBLE_FIGURE  Build a labeled preview composite from panel exports.
%
%   ebola_vhh.assemble_figure('fig5')
%   ebola_vhh.assemble_figure('fig5', 'RunPanels', true)

    run_panels = false;
    if ~isempty(varargin)
        for i = 1:2:numel(varargin)
            if strcmpi(varargin{i}, 'RunPanels')
                run_panels = logical(varargin{i+1});
            end
        end
    end

    M = ebola_vhh.load_manifest();
    F = M.figures.(fig_id);
    P = ebola_vhh.load_paths();
    panels = fieldnames(F.panels);

    if run_panels
        for i = 1:numel(panels)
            ebola_vhh.run_panel(fig_id, panels{i});
        end
    end

    layout = normalize_layout(F.layout);
    n_rows = numel(layout);
    n_cols = 1;
    for r = 1:n_rows
        n_cols = max(n_cols, numel(layout{r}));
    end

    fig = figure('Color', 'w', 'Name', F.id, ...
        'Position', [80 80 220*n_cols 200*n_rows]);
    tiled = tiledlayout(n_rows, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tiled, sprintf('%s.  %s', F.id, F.title), 'Interpreter', 'none');

    for r = 1:n_rows
        row = layout{r};
        if ~iscell(row)
            row = num2cell(row);
        end
        for c = 1:numel(row)
            pid = char(string(row{c}));
            nexttile((r-1)*n_cols + c);
            png = pick_panel_png(P.output, fig_id, pid);
            if isempty(png)
                axis off
                title(sprintf('%s (missing)', pid));
                continue
            end
            imshow(imread(png));
            title(pid, 'FontWeight', 'bold');
        end
        for c = (numel(row)+1):n_cols
            nexttile((r-1)*n_cols + c);
            axis off
        end
    end

    outdir = fullfile(P.output, fig_id, 'composites');
    if exist(outdir, 'dir') ~= 7
        mkdir(outdir);
    end
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, fullfile(outdir, [fig_id '_preview.png']), 'Resolution', 200);
        exportgraphics(fig, fullfile(outdir, [fig_id '_preview.pdf']), 'ContentType', 'image');
    else
        print(fig, fullfile(outdir, [fig_id '_preview.png']), '-dpng', '-r200');
    end
    fprintf('Wrote %s\n', fullfile(outdir, [fig_id '_preview.png']));
end

function rows = normalize_layout(layout)
    if iscell(layout) && ~isempty(layout) && iscell(layout{1})
        rows = layout;
        return
    end
    if iscell(layout)
        rows = cell(size(layout, 1), 1);
        for r = 1:size(layout, 1)
            rows{r} = layout(r, :);
        end
        return
    end
    rows = {layout};
end

function png = pick_panel_png(outroot, fig_id, pid)
    panel_dir = fullfile(outroot, fig_id, 'panels', pid);
    candidates = { ...
        fullfile(panel_dir, 'panel.png'), ...
        fullfile(panel_dir, 'cfd.png'), ...
        fullfile(panel_dir, 'vsv_intensity.png'), ...
        fullfile(panel_dir, 'vhh_per_virion.png')};
    png = '';
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file') == 2
            png = candidates{i};
            return
        end
    end
    d = dir(fullfile(panel_dir, '*.png'));
    if ~isempty(d)
        png = fullfile(d(1).folder, d(1).name);
    end
end
