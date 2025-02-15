clear;
load("RingSimu.mat");

%% PARAMETERS
saturationFactor = 0.6;
sz = 2;

color(1).hist = '#A1D6E4';
color(1).plot = '#2f8ca5';

color(2).hist = '#FFa67a';
color(2).plot = '#EF5000';

color(3).hist = '#CCFFCC';
color(3).plot = '#00a500';

%% DISPLAY
figure;

for i = 1:length(results)
    % excluding outliers
    locs = results(i).theta(:, ~results(i).outliers);
    errs = results(i).thetafinal(:, ~results(i).outliers);
    % locs = results(i).theta;
    % errs = results(i).thetafinal;
    
    % calculate pseudo color based on phi
    h = mod(rad2deg(locs(6,:)), 180) / 180;
    s = ones(size(h)) * saturationFactor;
    v = ones(size(h));
    orients = hsv2rgb([h(:), s(:), v(:)]);
    
    % calculate errors
    % error_xy = sqrt((errs(1, :).^2 + errs(2, :).^2) / 2);
    locs_r = sqrt((locs(1,:).^2 + locs(2,:).^2) / 2);
    object_r = sqrt((object(1, ~results(i).outliers).^2 + ...
        object(2, ~results(i).outliers).^2) / 2);
    error_r = locs_r - object_r;

    error_phi = rad2deg(errs(6,:));
    
    % fitting to normal distributions
    pd_r = fitdist(error_r', 'Normal');
    pd_phi = fitdist(error_phi', 'Normal');

    subplot(3, 3, i);
    scatter(locs(1,:), locs(2,:), sz, orients, 'filled');
    axis image;
    xlim([-100 100]); xticks([]);
    ylim([-100 100]); yticks([]);
    % axis off;
    set(gca, "Color", "#000000");

    % adding scale bar
    scalebar_length = 50;
    scalebar_x = 40;
    scalebar_y = -85;

    line([scalebar_x, scalebar_x + scalebar_length], [scalebar_y, scalebar_y], ...
        'Color', 'w', 'LineWidth', 2);

    subplot(3, 3, i+3);
    histogram(error_r, 'Normalization', 'pdf', 'BinWidth', 2, ...
        'FaceColor', color(i).hist);
    xlabel('Error_{xy}');
    title(['\sigma_{xy} = ' num2str(pd_r.sigma, '%.2f') ' nm']);
    hold on;
    xlim([-20 20]); xticks([-20 -10 0 10 20]);
    ylim([0 0.2]); yticks([0 0.1 0.2]);
    yr = normpdf(-20:20, pd_r.mu, pd_r.sigma);
    plot(-20:20, yr, 'LineWidth', 2, 'Color', color(i).plot);
    hold off;

    subplot(3, 3, i+6);
    histogram(error_phi, 'Normalization', 'pdf', 'BinWidth', 4, ...
        'FaceColor', color(i).hist);
    xlabel('Error_\phi');
    title(['\sigma_\phi = ' num2str(pd_phi.sigma, '%.2f') '°']);
    hold on;
    xlim([-40 40]); xticks([-40 -20 0 20 40]);
    ylim([0 .1]); yticks([0 .05 .1]);
    yr = normpdf(-40:40, pd_phi.mu, pd_phi.sigma);
    plot(-40:40, yr, 'LineWidth', 2, 'Color', color(i).plot);
    hold off;
end
