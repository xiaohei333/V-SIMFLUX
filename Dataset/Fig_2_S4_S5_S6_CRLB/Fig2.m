clear;
load("All_CRLB.mat");
load("Vortex_Simu_1.mat");
load("Vortex_Simu_2.mat");
load("V_SIMFLUX_Simu_1.mat");
load("V_SIMFLUX_Simu_2.mat");

g2_idx = 4;

%% EXTRACT
% check g2_idx
if mod(g2_idx, 1) || g2_idx < 1 || g2_idx > size(crlbv, 4)
    error('g2_idx must be integer and in the range [1,%d]', size(crlbv, 4));
end

crlbv_xy = squeeze(sqrt((crlbv(1,:,:,g2_idx).^2+crlbv(2,:,:,g2_idx).^2)/2));
crlbvs_xy = squeeze(sqrt((crlbvs(1,:,:,g2_idx).^2+crlbvs(2,:,:,g2_idx).^2)/2));

crlbv_phi = rad2deg(squeeze(crlbv(6,:,:,g2_idx)));
crlbvs_phi = rad2deg(squeeze(crlbvs(6,:,:,g2_idx)));

stdv = zeros(size(crlb_obj_v, [1 2 3]));
stdvs = stdv;
for i = 1:size(estv, 2)
    for j = 1:size(estv, 3)
        stdv(:,i,j) = get_statistics( ...
            squeeze(objv(:, i, j, g2_idx, :)), ...
            squeeze(estv(:, i, j, g2_idx, :)), ...
            squeeze(crlb_obj_v(:, i, j, g2_idx, :)), ...
            squeeze(outlierv(i, j, g2_idx, :))');
        stdvs(:,i,j) = get_statistics( ...
            squeeze(objvs(:, i, j, g2_idx, :)), ...
            squeeze(estvs(:, i, j, g2_idx, :)), ...
            squeeze(crlb_obj_vs(:, i, j, g2_idx, :)), ...
            squeeze(outliervs(i, j, g2_idx, :))');
    end
end

stdvxy = squeeze(sqrt((stdv(1,:,:).^2 + stdv(2,:,:).^2)/2));
stdvp = squeeze(stdv(6,:,:));
stdvsxy = squeeze(sqrt((stdvs(1,:,:).^2 + stdvs(2,:,:).^2)/2));
stdvsp = squeeze(stdvs(6,:,:));

lg = {'CRLB_{Vortex}', 'CRLB_{V-SIMFLUX}', '\sigma_{Vortex}', '\sigma_{V-SIMFLUX}'};

%% DISPLAY

figure;
subplot 221,
[h1, c1] = polarPcolor(theta, phi, crlbv_xy, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c1, 'CRLB_{xy} (nm)');
ax = gca;
set(ax, 'CLim', [1 8]);

subplot 222,
[h2, c2] = polarPcolor(theta, phi, crlbvs_xy, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c2, 'CRLB_{xy} (nm)');
ax = gca;
set(ax, 'CLim', [1 8]);

subplot 223,
[h3, c3] = polarPcolor(theta, phi, crlbv_phi, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c3, 'CRLB_{\phi} (°)');
ax = gca;
set(ax, 'CLim', [0 15]);

subplot 224,
[h4, c4] = polarPcolor(theta, phi, crlbvs_phi, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c4, 'CRLB_{\phi} (°)');
ax = gca;
set(ax, 'CLim', [0 15]);
sgtitle(['g_2 = ' num2str(g2(g2_idx), '%.2f')]);

figure,
subplot 221, plot(theta, mean(crlbv_xy, 1, 'omitmissing'), ...
    theta, mean(crlbvs_xy, 1, 'omitmissing'), 'LineWidth', 2);
hold on;
plot(20:10:90, mean(stdvxy, 1), ':ok', 'MarkerSize', 4);
plot(20:10:90, mean(stdvsxy, 1), ':^k', 'MarkerSize', 4);
hold off;
xlabel('Polar angle \theta (°)');
ylabel('CRLB_{xy} (nm)');
% legend(lg);
ylim([1 8]);
xlim([0 90]);
grid on;

subplot 222, plot(phi, mean(crlbv_xy(:, 21:end), 2, 'omitmissing'), ...
    phi, mean(crlbvs_xy(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
hold on;
plot(0:10:180, mean(stdvxy, 2), ':ok', 'MarkerSize', 4);
plot(0:10:180, mean(stdvsxy, 2), ':^k', 'MarkerSize', 4);
hold off;
xlabel('Azimuthal angle \phi (°)');
ylabel('CRLB_{xy} (nm)');
% legend(lg);
ylim([1 8]);
xlim([0 180]);
grid on;

subplot 223, plot(theta, mean(crlbv_phi, 1, 'omitmissing'), ...
    theta, mean(crlbvs_phi, 1, 'omitmissing'), 'LineWidth', 2);
hold on;
plot(20:10:90, rad2deg(mean(stdvp, 1)), ':ok', 'MarkerSize', 4);
plot(20:10:90, rad2deg(mean(stdvsp, 1)), ':^k', 'MarkerSize', 4);
hold off;
xlabel('Polar angle \theta (°)');
ylabel('CRLB_{\phi} (°)');
% legend(lg);
ylim([0 15]);
xlim([0 90]);
grid on;

subplot 224, plot(phi, mean(crlbv_phi(:, 21:10:end), 2, 'omitmissing'), ...
    phi, mean(crlbvs_phi(:, 21:10:end), 2, 'omitmissing'), 'LineWidth', 2);
hold on;
plot(0:10:180, rad2deg(mean(stdvp, 2)), ':ok', 'MarkerSize', 4);
plot(0:10:180, rad2deg(mean(stdvsp, 2)), ':^k', 'MarkerSize', 4);
hold off;
xlabel('Azimuthal angle \phi (°)');
ylabel('CRLB_{\phi} (°)');
legend(lg);
ylim([0 15]);
xlim([0 180]);
grid on;

sgtitle(['g_2 = ' num2str(g2(g2_idx), '%.2f')]);

%% PRINT
fprintf('XY enhancement: %.2f\n', ...
    mean(crlbv_xy(:, 21:end), 'all', 'omitmissing') / ...
    mean(crlbvs_xy(:, 21:end), 'all', 'omitmissing'));
fprintf('φ enhancement: %.2f\n', ...
    mean(crlbv_phi(:, 21:end), 'all', 'omitmissing') / ...
    mean(crlbvs_phi(:, 21:end), 'all', 'omitmissing'));
