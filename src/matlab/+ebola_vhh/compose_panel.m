function fig = compose_panel(fig_id, panel_id)
%COMPOSE_PANEL  Montage exported subpart PNGs into one panel preview.

    M = ebola_vhh.load_manifest();
    Pan = M.figures.(fig_id).panels.(panel_id);
    parts = ebola_vhh.as_cellstr(Pan.subparts);
    P = ebola_vhh.load_paths();
    panel_dir = fullfile(P.output, fig_id, 'panels', panel_id);

    imgs = {};
    labels = {};
    for i = 1:numel(parts)
        png = fullfile(panel_dir, [parts{i} '.png']);
        if exist(png, 'file') == 2
            imgs{end+1} = imread(png); %#ok<AGROW>
            labels{end+1} = parts{i}; %#ok<AGROW>
        end
    end

    if isempty(imgs)
        warning('ebola_vhh:noSubparts', ...
            'No exported PNGs for %s %s in %s', fig_id, panel_id, panel_dir);
        fig = [];
        return
    end

    fig = figure('Color', 'w', 'Name', sprintf('%s %s', fig_id, panel_id));
    n = numel(imgs);
    for i = 1:n
        ax = subplot(1, n, i);
        imshow(imgs{i}, 'Parent', ax);
        title(ax, labels{i}, 'Interpreter', 'none', 'FontSize', 9);
    end
    sgtitle(sprintf('%s%s  %s', fig_id, panel_id, Pan.title), ...
        'Interpreter', 'none', 'FontSize', 11);

    outdir = fullfile(P.output, fig_id, 'panels', panel_id);
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, fullfile(outdir, 'panel.png'), 'Resolution', 200);
        exportgraphics(fig, fullfile(outdir, 'panel.pdf'), 'ContentType', 'image');
    else
        print(fig, fullfile(outdir, 'panel.png'), '-dpng', '-r200');
    end
end
