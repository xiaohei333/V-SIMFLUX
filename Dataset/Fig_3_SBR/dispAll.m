load("SBR_Results.mat");
lw = 1.5;   % Line Width
sz = 4;     % Marker Size

%% EXTRACT
crlb_avg = squeeze(mean(mean(crlb, 2), 3));
est_std_avg = squeeze(mean(mean(est_std, 2), 3));
% est_bias_avg = squeeze(mean(mean(abs(est_bias), 2), 3));
% est_bias_avg = squeeze(mean(mean(est_bias, 2), 3));

crlb_avg(6:7, :) = rad2deg(crlb_avg(6:7, :));
est_std_avg(6:7, :) = rad2deg(est_std_avg(6:7, :));
% est_bias_avg(6:7, :) = rad2deg(est_bias_avg(6:7, :));

crlb_xy = sqrt((crlb_avg(1,:).^2 + crlb_avg(2,:).^2) / 2);
est_std_xy = sqrt((est_std_avg(1,:).^2 + est_std_avg(2,:).^2) / 2);
% est_bias_xy = sqrt((est_bias_avg(1,:).^2 + est_bias_avg(2,:).^2) / 2);

%%
figure;
subplot 131;
loglog(sbr, crlb_xy, '-', 'LineWidth', lw, 'Color', '#D95319');
hold on;
loglog(sbr, crlb_avg(3,:), '-', 'LineWidth', lw, 'Color', '#0072BD');
loglog(sbr, est_std_xy, ':ok', 'MarkerSize', sz);
loglog(sbr, est_std_avg(3,:), ':^k', 'MarkerSize', sz);
hold off;
xticks([10 1e2 1e3 1e4]); xlabel('SBR');
ylim([1e-2 1e3]); ylabel('Precision');
legend('CRLB_{xy} (nm)', 'CRLB_{z} (nm)', ...
    '\sigma_{xy} (nm)', '\sigma_{z} (nm)');
grid on;

subplot 132;
loglog(sbr, crlb_avg(6,:), '-', 'LineWidth', lw,'Color', '#77AC30');
hold on;
loglog(sbr, crlb_avg(7,:), '-', 'LineWidth', lw, 'Color', '#7E2F8E');
loglog(sbr, crlb_avg(8,:), '-', 'LineWidth', lw, 'Color', '#EDB120');
loglog(sbr, est_std_avg(6,:), ':ok', 'MarkerSize', sz);
loglog(sbr, est_std_avg(7,:), ':^k', 'MarkerSize', sz);
loglog(sbr, est_std_avg(8,:), ':>k', 'MarkerSize', sz);
hold off;
xticks([10 1e2 1e3 1e4]); xlabel('SBR');
ylim([1e-2 1e3]);
legend('CRLB_{\phi} (°)', 'CRLB_{\theta} (°)', 'CRLB_{g_2}', ...
    '\sigma_{\phi} (°)', '\sigma_{\theta} (°)', '\sigma_{g_2}');
grid on;

subplot 133;
loglog(sbr, crlb_avg(4,:), '-', 'LineWidth', lw, 'Color', '#4DBEEE');
hold on;
loglog(sbr, crlb_avg(5,:), '-', 'LineWidth', lw, 'Color', '#A2142F');
loglog(sbr, est_std_avg(4,:), ':ok', 'MarkerSize', sz);
loglog(sbr, est_std_avg(5,:), ':^k', 'MarkerSize', sz);
hold off;
xticks([10 1e2 1e3 1e4]); xlabel('SBR');
ylim([1e-2 1e3]);
legend('CRLB_{N}', 'CRLB_{bg}', ...
    '\sigma_{N}', '\sigma_{bg}');
grid on;

% figure;
% subplot 131;
% % semilogx(sbr, est_bias_xy./crlb_xy, sbr, est_bias_avg(3,:)./crlb_avg(3,:));
% semilogx(sbr, est_bias_avg(1,:)./crlb_avg(1,:), sbr, est_bias_avg(2,:)./crlb_avg(2,:), sbr, est_bias_avg(3,:)./crlb_avg(3,:));
% % semilogx(sbr, est_bias_avg(1,:), sbr, est_bias_avg(2,:), sbr, est_bias_avg(3,:));
%
% subplot 132;
% semilogx(sbr, est_bias_avg(6,:)./crlb_avg(6,:), sbr, est_bias_avg(7,:)./crlb_avg(7,:), sbr, est_bias_avg(8,:)./crlb_avg(8,:));
% % semilogx(sbr, est_bias_avg(6,:), sbr, est_bias_avg(7,:), sbr, est_bias_avg(8,:));
%
% subplot 133;
% semilogx(sbr, est_bias_avg(4,:)./crlb_avg(4,:), sbr, est_bias_avg(5,:)./crlb_avg(5,:));
% % semilogx(sbr, est_bias_avg(4,:), sbr, est_bias_avg(5,:));
