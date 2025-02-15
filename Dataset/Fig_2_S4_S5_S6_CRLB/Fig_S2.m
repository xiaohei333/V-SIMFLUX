clear;
load("All_CRLB.mat");

% dataName = 'V-SIMFLUX';
dataName = 'Vortex';

%% LIMITS
limFlag = true;    % set false to ignore limits below

xylim = [3 6];
zlim = [20 40];
thetalim = [2 8];
philim = [2 20];
g2lim = [0.02 0.1];

%% DISPLAY
switch lower(dataName)
    case 'vortex'
        crlb = crlbv;
    case 'v-simflux'
        crlb = crlbvs;
    otherwise
        error("The variable 'dataName' must be 'Vortex' or 'V-SIMFLUX'");
end

m = size(crlb, 4);
crlb_xy = squeeze(sqrt((crlb(1,:,:,:).^2 + crlbv(2,:,:,:).^2)/2));
crlb(6:7,:,:,:) = rad2deg(crlb(6:7,:,:,:));

figure;
for i = 1:m
    subplot(5, m, i);
    [h1(i), c1(i)] = polarPcolorNoDigits(theta, phi, crlb_xy(:,:,i), 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c1(i), 'CRLB_{xy} (nm)');
    % h1(i) = polarPcolorNoDigits(theta, phi, crlb_xy(:,:,i), 'Nspokes', 13, 'Ncircles', 4, 'colBar', 0);
    if limFlag
        ax = gca; set(ax, 'CLim', xylim);
    end

    subplot(5, m, i+m);
    [h2(i), c2(i)] = polarPcolorNoDigits(theta, phi, squeeze(crlb(3,:,:,i)), 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c2(i), 'CRLB_{z} (nm)');
    % h2(i) = polarPcolorNoDigits(theta, phi, squeeze(crlb(3,:,:,i)), 'Nspokes', 13, 'Ncircles', 4, 'colBar', 0);
    if limFlag
        ax = gca; set(ax, 'CLim', zlim);
    end

    subplot(5, m, i+m*2);
    [h3(i), c3(i)] = polarPcolorNoDigits(theta, phi, squeeze(crlb(7,:,:,i)), 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c3(i), 'CRLB_{\theta} (°)');
    % h3(i) = polarPcolorNoDigits(theta, phi, squeeze(crlb(7,:,:,i)), 'Nspokes', 13, 'Ncircles', 4, 'colBar', 0);
    if limFlag
        ax = gca; set(ax, 'CLim', thetalim);
    end

    subplot(5, m, i+m*3);
    [h4(i), c4(i)] = polarPcolorNoDigits(theta, phi, squeeze(crlb(6,:,:,i)), 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c4(i), 'CRLB_{\phi} (°)');
    % h4(i) = polarPcolorNoDigits(theta, phi, squeeze(crlb(6,:,:,i)), 'Nspokes', 13, 'Ncircles', 4, 'colBar', 0);
    if limFlag
        ax = gca; set(ax, 'CLim', philim);
    end

    subplot(5, m, i+m*4);
    [h5(i), c5(i)] = polarPcolorNoDigits(theta, phi, squeeze(crlb(8,:,:,i)), 'Nspokes', 13, 'Ncircles', 4);
    ylabel(c5(i), 'CRLB_{g_2}');
    % h5(i) = polarPcolorNoDigits(theta, phi, squeeze(crlb(8,:,:,i)), 'Nspokes', 13, 'Ncircles', 4, 'colBar', 0);
    if limFlag
        ax = gca; set(ax, 'CLim', g2lim);
    end
end
