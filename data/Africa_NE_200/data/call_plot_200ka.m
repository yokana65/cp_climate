

mean_val = mean(d(~isnan(d(:,2)), 2));
std_dev = std(d(~isnan(d(:,2)), 2));

% Define ylim based on 3-sigma range
y_min = mean_val - 5 * std_dev;
y_max = mean_val + 5 * std_dev;

a(nr) = axes('Position',[0.1 pos+.04 0.8 0.05],...
    'Visible','Off');
t(nr) = text(0.99,1.0,namestr,...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');

b(nr) = axes('Position',[0.1 pos 0.8 0.09],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'reverse',...
    "YAxisLocation","left", ...
    "ylim",[y_min, y_max]);





yt = d(:,2); % your time series y values
xt = d(:,1); % your time series x values
mean_yt = mean(d(~isnan(d(:,2)), 2));
hold on
fill_x = [xt, fliplr(xt)];
fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt<mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.1, 'EdgeAlpha', 0);
fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt>mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#0072BD', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);




line(d(:,1), d(:,2), "Color", "#000000", ...
    "LineWidth", 0.5, ...
    "Marker", ".", ...
    "MarkerSize", 4);
