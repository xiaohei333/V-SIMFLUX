clear;
load("Vortex_Simu_1.mat");
load("Vortex_Simu_2.mat");
load("V_SIMFLUX_Simu_1.mat");
load("V_SIMFLUX_Simu_2.mat");

theta = 20:10:90;
phi = 0:10:180;

lw = 2;
sz = 6;

lg = {'CRLB_{Vortex}', 'CRLB_{V-SIMFLUX}', '\sigma_{Vortex}', '\sigma_{V-SIMFLUX}'};

%% LIMITS
limFlag = true;    % set false to ignore limits below

xylim = [0 10];
zlim = [10 80];
thetalim = [0 10];
philim = [0 20];
g2lim = [0 0.1];

% xytick = [];

%% EXTRACT
crlbv = zeros(size(crlb_obj_v, [1 2 3]));
crlbvs = crlbv;
stdv = crlbv;
stdvs = crlbv;
for i = 1:size(estv, 2)
    for j = 1:size(estv, 3)
        for k = 1:size(estv, 4)
            [stdv(:,i,j,k), crlbv(:,i,j,k)] = get_statistics( ...
                squeeze(objv(:, i, j, k, :)), ...
                squeeze(estv(:, i, j, k, :)), ...
                squeeze(crlb_obj_v(:, i, j, k, :)), ...
                squeeze(outlierv(i, j, k, :))');
            [stdvs(:,i,j,k), crlbvs(:,i,j,k)] = get_statistics( ...
                squeeze(objvs(:, i, j, k, :)), ...
                squeeze(estvs(:, i, j, k, :)), ...
                squeeze(crlb_obj_vs(:, i, j, k, :)), ...
                squeeze(outliervs(i, j, k, :))');
        end
    end
end

%% DISPLAY
m = size(crlbv, 4);

crlbvxy = squeeze(sqrt((crlbv(1,:,:,:).^2 + crlbv(2,:,:,:).^2)/2));
crlbvsxy = squeeze(sqrt((crlbvs(1,:,:,:).^2 + crlbvs(2,:,:,:).^2)/2));
stdvxy = squeeze(sqrt((stdv(1,:,:,:).^2 + stdv(2,:,:,:).^2)/2));
stdvsxy = squeeze(sqrt((stdvs(1,:,:,:).^2 + stdvs(2,:,:,:).^2)/2));

crlbv(6:7,:,:,:) = rad2deg(crlbv(6:7,:,:,:));
crlbvs(6:7,:,:,:) = rad2deg(crlbvs(6:7,:,:,:));
stdv(6:7,:,:,:) = rad2deg(stdv(6:7,:,:,:));
stdvs(6:7,:,:,:) = rad2deg(stdvs(6:7,:,:,:));

figure;
for i = 1:m
    subplot(5, m, i);
    plot(theta, squeeze(mean(crlbvxy(:,:,i), 1)), ...
        theta, squeeze(mean(crlbvsxy(:,:,i), 1)), 'LineWidth', lw);
    hold on;
    plot(theta, squeeze(mean(stdvxy(:,:,i), 1)), ':ok', ...
        theta, squeeze(mean(stdvsxy(:,:,i), 1)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([20 90]); xticks([20 55 90]); grid on;
    if i == 1
        ylabel('Precsion_{xy} (nm)');
    end
    if limFlag
        ylim(xylim);
        % yticks(xytick);
    end

    subplot(5, m, i+m);
    plot(theta, squeeze(mean(crlbv(3,:,:,i), 2)), ...
        theta, squeeze(mean(crlbvs(3,:,:,i), 2)), 'LineWidth', lw);
    hold on;
    plot(theta, squeeze(mean(stdv(3,:,:,i), 2)), ':ok', ...
        theta, squeeze(mean(stdvs(3,:,:,i), 2)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([20 90]); xticks([20 55 90]); grid on;
    if i == 1
        ylabel('Precsion_{z} (nm)');
    end
    if limFlag
        ylim(zlim);
        % yticks(ztick);
    end

    subplot(5, m, i+m*2);
    plot(theta, squeeze(mean(crlbv(6,:,:,i), 2)), ...
        theta, squeeze(mean(crlbvs(6,:,:,i), 2)), 'LineWidth', lw);
    hold on;
    plot(theta, squeeze(mean(stdv(6,:,:,i), 2)), ':ok', ...
        theta, squeeze(mean(stdvs(6,:,:,i), 2)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([20 90]); xticks([20 55 90]); grid on;
    if i == 1
        ylabel('Precsion_{\phi} (°)');
    end
    if limFlag
        ylim(philim);
        % yticks(phitick);
    end

    subplot(5, m, i+m*3);
    plot(theta, squeeze(mean(crlbv(7,:,:,i), 2)), ...
        theta, squeeze(mean(crlbvs(7,:,:,i), 2)), 'LineWidth', lw);
    hold on;
    plot(theta, squeeze(mean(stdv(7,:,:,i), 2)), ':ok', ...
        theta, squeeze(mean(stdvs(7,:,:,i), 2)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([20 90]); xticks([20 55 90]); grid on;
    if i == 1
        ylabel('Precsion_{\theta} (°)');
    end
    if limFlag
        ylim(thetalim);
        % yticks(thetatick);
    end

    subplot(5, m, i+m*4);
    plot(theta, squeeze(mean(crlbv(8,:,:,i), 2)), ...
        theta, squeeze(mean(crlbvs(8,:,:,i), 2)), 'LineWidth', lw);
    hold on;
    plot(theta, squeeze(mean(stdv(8,:,:,i), 2)), ':ok', ...
        theta, squeeze(mean(stdvs(8,:,:,i), 2)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([20 90]); xticks([20 55 90]); grid on;
    xlabel('Polar Angle \theta (°)');
    if i == 1
        ylabel('Precsion_{g2}');
    end
    if limFlag
        ylim(g2lim);
        % yticks(g2tick);
    end
    
end

figure;
for i = 1:m
    subplot(5, m, i);
    plot(phi, squeeze(mean(crlbvxy(:,:,i), 2)), ...
        phi, squeeze(mean(crlbvsxy(:,:,i), 2)), 'LineWidth', lw);
    hold on;
    plot(phi, squeeze(mean(stdvxy(:,:,i), 2)), ':ok', ...
        phi, squeeze(mean(stdvsxy(:,:,i), 2)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([0 180]); xticks([0 45 90 135 180]); grid on;
    if i == 1
        ylabel('Precsion_{xy} (nm)');
    end
    if limFlag
        ylim(xylim);
        % yticks(xytick);
    end

    subplot(5, m, i+m);
    plot(phi, squeeze(mean(crlbv(3,:,:,i), 3)), ...
        phi, squeeze(mean(crlbvs(3,:,:,i), 3)), 'LineWidth', lw);
    hold on;
    plot(phi, squeeze(mean(stdv(3,:,:,i), 3)), ':ok', ...
        phi, squeeze(mean(stdvs(3,:,:,i), 3)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([0 180]); xticks([0 45 90 135 180]); grid on;
    if i == 1
        ylabel('Precsion_{z} (nm)');
    end
    if limFlag
        ylim(zlim);
        % yticks(ztick);
    end

    subplot(5, m, i+m*2);
    plot(phi, squeeze(mean(crlbv(6,:,:,i), 3)), ...
        phi, squeeze(mean(crlbvs(6,:,:,i), 3)), 'LineWidth', lw);
    hold on;
    plot(phi, squeeze(mean(stdv(6,:,:,i), 3)), ':ok', ...
        phi, squeeze(mean(stdvs(6,:,:,i), 3)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([0 180]); xticks([0 45 90 135 180]); grid on;
    if i == 1
        ylabel('Precsion_{\phi} (°)');
    end
    if limFlag
        ylim(philim);
        % yticks(phitick);
    end

    subplot(5, m, i+m*3);
    plot(phi, squeeze(mean(crlbv(7,:,:,i), 3)), ...
        phi, squeeze(mean(crlbvs(7,:,:,i), 3)), 'LineWidth', lw);
    hold on;
    plot(phi, squeeze(mean(stdv(7,:,:,i), 3)), ':ok', ...
        phi, squeeze(mean(stdvs(7,:,:,i), 3)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([0 180]); xticks([0 45 90 135 180]); grid on;
    if i == 1
        ylabel('Precsion_{\theta} (°)');
    end
    if limFlag
        ylim(thetalim);
        % yticks(thetatick);
    end

    subplot(5, m, i+m*4);
    plot(phi, squeeze(mean(crlbv(8,:,:,i), 3)), ...
        phi, squeeze(mean(crlbvs(8,:,:,i), 3)), 'LineWidth', lw);
    hold on;
    plot(phi, squeeze(mean(stdv(8,:,:,i), 3)), ':ok', ...
        phi, squeeze(mean(stdvs(8,:,:,i), 3)), ':^k', 'MarkerSize', sz);
    hold off;
    xlim([0 180]); xticks([0 45 90 135 180]); grid on;
    xlabel('Azimuthal Angle \phi (°)');
    if i == 1
        ylabel('Precsion_{g2}');
    end
    if limFlag
        ylim(g2lim);
        % yticks(g2tick);
    end
    
end
