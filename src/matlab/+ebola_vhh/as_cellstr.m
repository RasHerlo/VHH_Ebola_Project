function parts = as_cellstr(x)
%AS_CELLSTR  Normalize JSON-decoded string lists (1-element arrays become char).

    if ischar(x)
        parts = {x};
    elseif isstring(x)
        parts = cellstr(x(:));
    elseif iscell(x)
        parts = cellfun(@char, x(:), 'UniformOutput', false);
    else
        error('ebola_vhh:as_cellstr', 'Cannot convert type %s', class(x));
    end
end
