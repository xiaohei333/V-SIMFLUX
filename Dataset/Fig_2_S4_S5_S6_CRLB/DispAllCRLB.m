clear;
load("All_CRLB.mat");

%% OPTIONS
limFlag = false;    % set false to ignore limits below
g2_idx = 4;

xylim = [1 8];
zlim = [20 40];
thetalim = [2 4];
philim = [0 15];
g2lim = [0.02 0.1];
% Nlim = [100 250];
% bglim = [0.28 0.36];

xytick = [2 4 6 8];
ztick = [20 30 40];
phitick = [0 5 10 15];
thetatick = [2 3 4];
g2tick = [0.02 0.04 0.06 0.08 0.1];

lg = {'Vortex', 'V-SIMFLUX'};

%% XY
% check g2_idx
if mod(g2_idx, 1) || g2_idx < 1 || g2_idx > size(crlbv, 4)
    error('g2_idx must be integer and in the range [1,%d]', size(crlbv, 4));
end

crlbv_xy = squeeze(sqrt((crlbv(1,:,:,g2_idx).^2 + crlbv(2,:,:,g2_idx).^2)/2));
crlbvs_xy = squeeze(sqrt((crlbvs(1,:,:,g2_idx).^2 + crlbvs(2,:,:,g2_idx).^2)/2));
% crbv_xy = squeeze((crbv(1,:,:,g2_idx) + crbv(2,:,:,g2_idx))/2);
% crbvs_xy = squeeze((crbvs(1,:,:,g2_idx) + crbvs(2,:,:,g2_idx))/2);

figure;
subplot 221,
[h1, c1] = polarPcolor(theta, phi, crlbv_xy, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c1, 'CRLB_{xy} (nm)');
if limFlag
    ax = gca; set(ax, 'CLim', xylim);
end

subplot 222,
[h2, c2] = polarPcolor(theta, phi, crlbvs_xy, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c2, 'CRLB_{xy} (nm)');
if limFlag
    ax = gca; set(ax, 'CLim', xylim);
end

subplot 223,
plot(theta, mean(crlbv_xy, 1, 'omitmissing'), ...
    theta, mean(crlbvs_xy, 1, 'omitmissing'), 'LineWidth', 2);
xlabel('Polar angle \theta (°)');
ylabel('CRLB_{xy} (nm)');
% legend(lg);
xlim([0 90]); xticks([0 30 60 90]);
if limFlag
    ylim(xylim); yticks(xytick);
end

subplot 224,
plot(phi, mean(crlbv_xy(:, 21:end), 2, 'omitmissing'), ...
    phi, mean(crlbvs_xy(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
xlabel('Azimuthal angle \phi (°)');
ylabel('CRLB_{xy} (nm)');
legend(lg);
xlim([0 180]); xticks([0 45 90 135 180]);
if limFlag
    ylim(xylim); yticks(xytick);
end

xyfold = mean(crlbv_xy(:, 21:end), 'all', 'omitmissing') / ...
    mean(crlbvs_xy(:, 21:end), 'all', 'omitmissing');
sgtitle(['XY enhance: ' sprintf('%.2f', xyfold) '-fold (\theta\geq20°) ' ...
    'with g_2 = ' sprintf('%.2f', g2(g2_idx))]);

%% Z
crbv_z = squeeze(crlbv(3,:,:,g2_idx));
crbvs_z = squeeze(crlbvs(3,:,:,g2_idx));

figure,
subplot 221,
[~,c1] = polarPcolor(theta, phi, crbv_z, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c1, 'CRLB_{z} (nm)');
if limFlag
    ax = gca; set(ax, 'CLim', zlim);
end

subplot 222,
[~, c2] = polarPcolor(theta, phi, crbvs_z, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c2, 'CRLB_{z} (nm)');
if limFlag
    ax = gca; set(ax, 'CLim', zlim);
end

subplot 223,
plot(theta, mean(crbv_z, 1, 'omitmissing'), ...
    theta, mean(crbvs_z, 1, 'omitmissing'), 'LineWidth', 2);
xlabel('Polar angle \theta (°)');
ylabel('CRLB_{z} (nm)');
% legend(lg);
xlim([0 90]); xticks([0 30 60 90]);
if limFlag
    ylim(zlim); yticks(ztick);
end

subplot 224,
plot(phi, mean(crbv_z(:, 21:end), 2, 'omitmissing'), ...
    phi, mean(crbvs_z(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
xlabel('Azimuthal angle \phi (°)');
ylabel('CRLB_{z} (nm)');
% legend(lg);
xlim([0 180]); xticks([0 45 90 135 180]);
if limFlag
    ylim(zlim); yticks(ztick);
end

zfold = mean(crbv_z(:, 21:end), 'all', 'omitmissing') / ...
    mean(crbvs_z(:, 21:end), 'all', 'omitmissing');
sgtitle(['Z enhance: ' sprintf('%.2f', zfold) '-fold (\theta\geq20°) ' ...
    'with g_2 = ' sprintf('%.2f', g2(g2_idx))]);

%% THETA
crbv_theta = rad2deg(squeeze(crlbv(7,:,:,g2_idx)));
crbvs_theta = rad2deg(squeeze(crlbvs(7,:,:,g2_idx)));

figure,
subplot 221,
[~,c1] = polarPcolor(theta, phi, crbv_theta, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c1, 'CRLB_{\theta} (°)');
if limFlag
    ax = gca; set(ax, 'CLim', thetalim);
end

subplot 222,
[~, c2] = polarPcolor(theta, phi, crbvs_theta, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c2, 'CRLB_{\theta} (°)');
if limFlag
    ax = gca; set(ax, 'CLim', thetalim);
end

subplot 223,
plot(theta, mean(crbv_theta, 1, 'omitmissing'), ...
    theta, mean(crbvs_theta, 1, 'omitmissing'), 'LineWidth', 2);
xlabel('Polar angle \theta (°)');
ylabel('CRLB_{\theta} (°)');
% legend(lg);
xlim([0 90]); xticks([0 30 60 90]);
if limFlag
    ylim(thetalim); yticks(thetatick);
end

subplot 224,
plot(phi, mean(crbv_theta(:, 21:end), 2, 'omitmissing'), ...
    phi, mean(crbvs_theta(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
xlabel('Azimuthal angle \phi (°)');
ylabel('CRLB_{\theta} (°)');
% legend(lg);
xlim([0 180]); xticks([0 45 90 135 180]);
if limFlag
    ylim(thetalim); yticks(thetatick);
end

thetafold = mean(crbv_theta(:, 21:end), 'all', 'omitmissing') / ...
    mean(crbvs_theta(:, 21:end), 'all', 'omitmissing');
sgtitle(['\theta enhance: ' sprintf('%.2f', thetafold) '-fold (\theta\geq20°) ' ...
    'with g_2 = ' sprintf('%.2f', g2(g2_idx))]);

%% PHI
crbv_phi = rad2deg(squeeze(crlbv(6,:,:,g2_idx)));
crbvs_phi = rad2deg(squeeze(crlbvs(6,:,:,g2_idx)));

figure,
subplot 221,
[~,c1] = polarPcolor(theta, phi, crbv_phi, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c1, 'CRLB_{\phi} (°)');
if limFlag
    ax = gca; set(ax, 'CLim', philim);
end

subplot 222,
[~, c2] = polarPcolor(theta, phi, crbvs_phi, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c2, 'CRLB_{\phi} (°)');
if limFlag
    ax = gca; set(ax, 'CLim', philim);
end

subplot 223,
plot(theta, mean(crbv_phi, 1, 'omitmissing'), ...
    theta, mean(crbvs_phi, 1, 'omitmissing'), 'LineWidth', 2);
xlabel('Polar angle \theta (°)'); ylabel('CRLB_{\phi} (°)');
% legend(lg);
xlim([0 90]); xticks([0 30 60 90]);
if limFlag
    ylim(philim); yticks(phitick);
end

subplot 224,
plot(phi, mean(crbv_phi(:, 21:end), 2, 'omitmissing'), ...
    phi, mean(crbvs_phi(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
xlabel('Azimuthal angle \phi (°)'); ylabel('CRLB_{\phi} (°)');
% legend(lg);
xlim([0 180]); xticks([0 45 90 135 180]);
if limFlag
    ylim(philim); yticks(phitick);
end

phifold = mean(crbv_phi(:, 21:end), 'all', 'omitmissing') / ...
    mean(crbvs_phi(:, 21:end), 'all', 'omitmissing');
sgtitle(['\phi enhance: ' sprintf('%.2f', phifold) '-fold (\theta\geq20°) ' ...
    'with g_2 = ' sprintf('%.2f', g2(g2_idx))]);

%% G2
crbv_g = squeeze(crlbv(8,:,:,g2_idx));
crbvs_g = squeeze(crlbvs(8,:,:,g2_idx));

figure,
subplot 221,
[~,c1] = polarPcolor(theta, phi, crbv_g, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c1, 'CRLB_{g_2}');
if limFlag
    ax = gca; set(ax, 'CLim', g2lim);
end

subplot 222,
[~, c2] = polarPcolor(theta, phi, crbvs_g, 'Nspokes', 13, 'Ncircles', 4);
ylabel(c2, 'CRLB_{g_2}');
if limFlag
    ax = gca; set(ax, 'CLim', g2lim);
end

subplot 223,
plot(theta, mean(crbv_g, 1, 'omitmissing'), ...
    theta, mean(crbvs_g, 1, 'omitmissing'), 'LineWidth', 2);
xlabel('Polar angle \theta (°)'); ylabel('CRLB_{g_2}');
% legend(lg);
xlim([0 90]); xticks([0 30 60 90]);
if limFlag
    ylim(g2lim); yticks(g2tick);
end

subplot 224,
plot(phi, mean(crbv_g(:, 21:end), 2, 'omitmissing'), ...
    phi, mean(crbvs_g(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
xlabel('Azimuthal angle \phi (°)'); ylabel('CRLB_{g_2}');
% legend(lg);
xlim([0 180]); xticks([0 45 90 135 180]);
if limFlag
    ylim(g2lim); yticks(g2tick);
end

g2fold = mean(crbv_g(:, 21:end), 'all', 'omitmissing') / ...
    mean(crbvs_g(:, 21:end), 'all', 'omitmissing');
sgtitle(['g_2 enhance: ' sprintf('%.2f', g2fold) '-fold (\theta\geq20°) ' ...
    'with g_2 = ' sprintf('%.2f', g2(g2_idx))]);

%% N
if ~limFlag || exist('Nlim', 'var')
    crbv_N = squeeze(crlbv(4,:,:,g2_idx));
    crbvs_N = squeeze(crlbvs(4,:,:,g2_idx));

    figure,
    subplot 221,
    [~,c1] = polarPcolor(theta, phi, crbv_N, 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c1, 'CRLB_{N}');
    if limFlag
        ax = gca; set(ax, 'CLim', Nlim);
    end

    subplot 222,
    [~, c2] = polarPcolor(theta, phi, crbvs_N, 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c2, 'CRLB_{N}');
    if limFlag
        ax = gca; set(ax, 'CLim', Nlim);
    end

    subplot 223,
    plot(theta, mean(crbv_N, 1, 'omitmissing'), ...
        theta, mean(crbvs_N, 1, 'omitmissing'), 'LineWidth', 2);
    xlabel('Polar angle \theta (°)'); ylabel('CRLB_{N}');
    % legend(lg);
    xlim([0 90]); xticks([0 30 60 90]);
    if limFlag
        ylim(Nlim); yticks(Ntick);
    end

    subplot 224,
    plot(phi, mean(crbv_N(:, 21:end), 2, 'omitmissing'), ...
        phi, mean(crbvs_N(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
    xlabel('Azimuthal angle \phi (°)'); ylabel('CRLB_{N}');
    % legend(lg);
    xlim([0 180]); xticks([0 45 90 135 180]);
    if limFlag
        ylim(Nlim); yticks(Ntick);
    end

    Nfold = mean(crbv_N(:, 21:end), 'all', 'omitmissing') / ...
        mean(crbvs_N(:, 21:end), 'all', 'omitmissing');
    sgtitle(['N enhance: ' sprintf('%.2f', Nfold) '-fold (\theta\geq20°) ' ...
    'with g_2 = ' sprintf('%.2f', g2(g2_idx))]);
end

%% BG
if ~limFlag || exist('Nlim', 'var')
    crbv_bg = squeeze(crlbv(5,:,:,g2_idx));
    crbvs_bg = squeeze(crlbvs(5,:,:,g2_idx));

    figure,
    subplot 221,
    [~,c1] = polarPcolor(theta, phi, crbv_bg, 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c1, 'CRLB_{bg}');
    if limFlag
        ax = gca; set(ax, 'CLim', bglim);
    end

    subplot 222,
    [~, c2] = polarPcolor(theta, phi, crbvs_bg, 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c2, 'CRLB_{bg}');
    if limFlag
        ax = gca; set(ax, 'CLim', bglim);
    end

    subplot 223,
    plot(theta, mean(crbv_bg, 1, 'omitmissing'), ...
        theta, mean(crbvs_bg, 1, 'omitmissing'), 'LineWidth', 2);
    xlabel('Polar angle \theta (°)'); ylabel('CRLB_{bg}');
    % legend(lg);
    xlim([0 90]); xticks([0 30 60 90]);
    if limFlag
        ylim(bglim); yticks(bgtick);
    end

    subplot 224,
    plot(phi, mean(crbv_bg(:, 21:end), 2, 'omitmissing'), ...
        phi, mean(crbvs_bg(:, 21:end), 2, 'omitmissing'), 'LineWidth', 2);
    xlabel('Azimuthal angle \phi (°)'); ylabel('CRLB_{bg}');
    % legend(lg);
    xlim([0 90]); xticks([0 45 90 135 180]);
    if limFlag
        ylim(bglim); yticks(bgtick);
    end

    bgfold = mean(crbv_bg(:, 21:end), 'all', 'omitmissing') / ...
        mean(crbvs_bg(:, 21:end), 'all', 'omitmissing');
    sgtitle(['BG enhance: ' sprintf('%.2f', bgfold) '-fold (\theta\geq20°) ' ...
    'with g_2 = ' sprintf('%.2f', g2(g2_idx))]);
end
