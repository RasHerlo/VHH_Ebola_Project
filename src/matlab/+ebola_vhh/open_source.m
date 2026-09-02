function fig = open_source(src)
%OPEN_SOURCE  Open a legacy .fig or raster image as a MATLAB figure.

    if exist(src, 'file') ~= 2
        error('ebola_vhh:missingSource', 'Source file not found:\n  %s', src);
    end

    [~, ~, ext] = fileparts(src);
    switch lower(ext)
        case '.fig'
            fig = openfig(src, 'new', 'visible');
        case {'.tif', '.tiff', '.jpg', '.jpeg', '.png', '.bmp'}
            fig = figure('Color', 'w');
            ax = axes(fig);
            img = imread(src);
            imshow(img, 'Parent', ax);
            axis(ax, 'image');
        otherwise
            error('ebola_vhh:unsupportedSource', ...
                'Cannot open %s as a MATLAB figure.', src);
    end
end
