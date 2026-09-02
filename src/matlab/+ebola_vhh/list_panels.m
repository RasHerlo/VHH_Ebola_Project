function list_panels()
%LIST_PANELS  Print the figure/panel tree from the manifest.

    M = ebola_vhh.load_manifest();
    figs = fieldnames(M.figures);
    fprintf('\n%s\n', M.scope);
    fprintf('%s\n\n', M.numbering_note);

    for i = 1:numel(figs)
        F = M.figures.(figs{i});
        fprintf('== %s  %s ==\n', F.id, F.title);
        fprintf('   %s\n', F.story);
        fprintf('   assemble: ebola_vhh.%s.assemble\n', F.id);
        panels = fieldnames(F.panels);
        for p = 1:numel(panels)
            Pan = F.panels.(panels{p});
            fprintf('   [%s] %s\n', panels{p}, Pan.title);
            fprintf('       %s\n', Pan.claim);
            fprintf('       %s  subparts: %s\n', Pan.call, strjoin(ebola_vhh.as_cellstr(Pan.subparts), ', '));
        end
        fprintf('\n');
    end

    if isfield(M, 'candidates')
        fprintf('== candidates (not in the draft) ==\n');
        cands = fieldnames(M.candidates);
        for i = 1:numel(cands)
            C = M.candidates.(cands{i});
            fprintf('   %s  %s\n       %s\n', cands{i}, C.call, C.claim);
        end
        fprintf('\n');
    end
end
