clear;
load("results.mat");

sz = 2.5;

color(1).plot = "#1a6fe4";
color(1).hist = "#7caffa";
color(2).plot = "#f5403d";
color(2).hist = "#f77276";

%% DISPLAY
results(6:7,:,:) = rad2deg(results(6:7,:,:));

figure;
for i = 1:size(results, 3)
    % fit to normal distributions
    pd_x(i) = fitdist(results(1,:,i)', 'Normal');
    pd_y(i) = fitdist(results(2,:,i)', 'Normal');
    pd_p(i) = fitdist(results(6,:,i)', 'Normal');
    pd_t(i) = fitdist(results(7,:,i)', 'Normal');

    subplot(4, 2, [1 3]);
    hold on;
    plot(results(1,:,i), results(2,:,i), 'o', ...
        'MarkerEdgeColor', color(i).plot, ...
        'MarkerFaceColor', color(i).plot, ...
        'MarkerSize', sz);
    axis square; box on;
    hold off;

    subplot(4, 2, [2 4]);
    hold on;
    plot(results(6,:,i), results(7,:,i), 'o', ...
        'MarkerEdgeColor', color(i).plot, ...
        'MarkerFaceColor', color(i).plot, ...
        'MarkerSize', sz);
    axis square; box on;
    hold off;

    x = -20:0.1:20;
    yx = normpdf(x, pd_x(i).mu, pd_x(i).sigma);
    yy = normpdf(x, pd_y(i).mu, pd_y(i).sigma);

    subplot(4, 2, 5);
    hold on;
    plot(x, yx, 'Color', color(i).hist, 'LineWidth', 1.5);
    fill([x, fliplr(x)], [yx, zeros(size(yx))], 'b', 'FaceColor', color(i).hist, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold off;


    subplot(4, 2, 7);
    hold on;
    plot(x, yy, 'Color', color(i).hist, 'LineWidth', 1.5);
    fill([x, fliplr(x)], [yy, zeros(size(yy))], 'b', 'FaceColor', color(i).hist, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xlim([-20 20]); ylim([0 0.3]);
    hold off;
    xlabel('y (nm)');
    ylabel('Probability density');

    p = 30:0.1:60;
    yp = normpdf(p, pd_p(i).mu, pd_p(i).sigma);
    t = 30:0.1:60;
    yt = normpdf(t, pd_t(i).mu, pd_t(i).sigma);

    subplot(4, 2, 6);
    hold on;
    plot(p, yp, 'Color', color(i).hist, 'LineWidth', 1.5);
    fill([p, fliplr(p)], [yp, zeros(size(yp))], 'b', 'FaceColor', color(i).hist, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold off;


    subplot(4, 2, 8);
    hold on;
    plot(t, yt, 'Color', color(i).hist, 'LineWidth', 1.5);
    fill([t, fliplr(t)], [yt, zeros(size(yt))], 'b', 'FaceColor', color(i).hist, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold off;

end

subplot(4, 2, [1 3]);
xlim([-20 20]); xticks([-20 0 20]);
ylim([-20 20]); yticks([-20 0 20]);
xlabel('x (nm)');
ylabel('y (nm)');

subplot(4, 2, [2 4]);
xlim([30 60]); xticks([30 45 60]);
ylim([30 60]); yticks([30 45 60]);
xlabel('\phi (°)');
ylabel('\theta (°)');

subplot(4, 2, 5);
xlim([-20 20]); ylim([0 0.3]);
xlabel('x (nm)');
ylabel('Probability density');
title(['X: ' num2str(pd_x(1).sigma/pd_x(2).sigma, '%.2f') '-fold']);

subplot(4, 2, 7);
xlim([-20 20]); ylim([0 0.3]);
xlabel('y (nm)');
ylabel('Probability density');
title(['Y: ' num2str(pd_y(1).sigma/pd_y(2).sigma, '%.2f') '-fold']);

subplot(4, 2, 6);
xlim([30 60]); ylim([0 0.3]);
xlabel('\phi (°)');
ylabel('Probability density');
title(['\phi: ' num2str(pd_p(1).sigma/pd_p(2).sigma, '%.2f') '-fold']);

subplot(4, 2, 8);
xlim([30 60]); ylim([0 0.3]);
xlabel('\theta (°)');
ylabel('Probability density');
title(['\theta: ' num2str(pd_t(1).sigma/pd_t(2).sigma, '%.2f') '-fold']);
