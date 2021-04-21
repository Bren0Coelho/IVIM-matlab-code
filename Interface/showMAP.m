function [im, c] = showMAP(parametric_map, label)
    im  = image(parametric_map, 'CDataMapping', 'scaled');
    c = colorbar;
    c.Label.String = label;
    c.Label.FontSize = 12;
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'FontSize', 12);
end