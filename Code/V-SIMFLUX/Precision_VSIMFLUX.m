close all;
clearvars;

% set parameters
params = set_parameters_vortex_sim;
params.Ncfg = 1000; % generate and fit Ncfg randomized model PSFs.
params.flg_parallel = 1;

[wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);

SIM.L = 2; % 条纹方向
SIM.K = 3; % 相移步数
SIM.Betal = [pi/2 pi]; % a global angular offset
SIM.Chi = [0 0]; % size(Chi) = 1*L
SIM.lp = 561/2/1.49; % 条纹周期（nm）
SIM.m = [0.9 0.9]; % 调制深度

object = zeros(params.numparams,params.Ncfg);
allpsfs = zeros(params.Mx,params.My,params.K,SIM.L,SIM.K,params.Ncfg);
dallpsfsdtheta = zeros(params.Mx,params.My,params.K,params.numparams,SIM.L,SIM.K,params.Ncfg);

for jj = 1:params.Ncfg

    % true parameters
    dx = 0;
    dy = 0;
    dz = 0;
    Nphotons = 4000;
    Nbackground = 10;
    dazim = pi/4;
    dpola = pi/2;
    dg2 = 0.75;

    ROIxy = [0 0];

    object(:,jj) = [dx dy dz Nphotons Nbackground dazim dpola dg2];

    [allpsfs(:,:,params.K,:,:,jj) , dallpsfsdtheta(:,:,:,:,:,:,jj)] ...
  = poissonrate_VSIMFLUX(params,SIM,object(:,jj),PupilMatrix,allzernikes,wavevector,wavevectorzimm,ROIxy);

end


allspots = 1e12*imnoise(allpsfs*1e-12,'poisson');



figure,
for l = 1:2
    for k = 1:3
        subplot(2,3,(l-1)*3+k);
        imagesc(squeeze(allpsfs(:,:,1,l,k)));
        axis image;
    end
end
% figure,
% imagesc(squeeze(sum(allpsfs,[4 5])));
% axis image;

%% MLE fitting routine
allpsfs_V(:,:,1,:) = squeeze(sum(allspots,[4 5]));
[thetainit_V] = initialvalues(allpsfs_V,params);
[thetastore_V,~,~,~,~] = localization(allpsfs_V,thetainit_V,params);
theta_V = thetastore_V(:,:,end);
thetainit_VSIMFLUX = theta_V;

roixy0 = zeros(params.Ncfg,2);
[thetastore_VSIMFLUX,mu,dmudtheta,merit,numiters] = localization_VSIMFLUX(allspots,thetainit_VSIMFLUX,params,SIM,roixy0);
theta_VSIMFLUX = thetastore_VSIMFLUX(:,:,end);


figure;
plot(theta_V(1,:),theta_V(2,:),'.'); hold on;
plot(theta_VSIMFLUX(1,:),theta_VSIMFLUX(2,:),'.'); 
axis image;

tmpcfg_V = theta_V(end-2,:)>pi;
theta_V(end-2,tmpcfg_V) = theta_V(end-2,tmpcfg_V)-pi;
theta_V(end-1,tmpcfg_V) = pi-theta_V(end-1,tmpcfg_V);

tmpcfg_VSIMFLUX = theta_VSIMFLUX(end-2,:)>pi;
theta_VSIMFLUX(end-2,tmpcfg_VSIMFLUX) = theta_VSIMFLUX(end-2,tmpcfg_VSIMFLUX)-pi;
theta_VSIMFLUX(end-1,tmpcfg_VSIMFLUX) = pi-theta_VSIMFLUX(end-1,tmpcfg_VSIMFLUX);

% figure;
% subplot(1,2,1);
% plot(theta_V(6,:)); hold on;
% plot(theta_V(7,:)); ylim([0 pi]);
% subplot(1,2,2);
% plot(theta_VSIMFLUX(6,:)); hold on;
% plot(theta_VSIMFLUX(7,:)); ylim([0 pi]);

% figure;
% subplot(1,2,1);
% histogram(theta_V(6,:)); hold on;
% histogram(theta_V(7,:)); 
% subplot(1,2,2);
% histogram(theta_VSIMFLUX(6,:)); hold on;
% histogram(theta_VSIMFLUX(7,:));

figure;
subplot(1,2,1);
histogram(theta_V(6,:)); hold on;
histogram(theta_VSIMFLUX(6,:)); 
subplot(1,2,2);
histogram(theta_V(7,:)); hold on;
histogram(theta_VSIMFLUX(7,:));


%% Data analysis

allpsfs_Vortex(:,:,1,:) = squeeze(sum(allpsfs,[4 5]));
dallpsfsdtheta_Vortex(:,:,1,:,:) = squeeze(sum(dallpsfsdtheta,[5 6]));
[crlb_Vortex,~] = get_fisher_crlb(params,allpsfs_Vortex,dallpsfsdtheta_Vortex);
mean(crlb_Vortex(3,:))
std(theta_V(3,:))

% [crlb_VSIMFLUX,rcondstore] = get_fisher_crlb_VSIMFLUX(params,SIM,allpsfs,dallpsfsdtheta);
[crlb_VSIMFLUX,rcondstore] = get_fisher_crlb_VSIMFLUX(params,SIM,mu,dmudtheta);
[outliers] = get_outliers(theta_VSIMFLUX,merit,numiters,params);
[thetafinal,thetamean,thetastd,crlbmean] = get_statistics(params,object,theta_VSIMFLUX,crlb_VSIMFLUX,outliers);

% figure,
% plot(thetafinal(1,:),thetafinal(2,:),'.'); 
% axis image;
% 
% figure,
% plot(thetafinal(6,:)); hold on;
% plot(thetafinal(7,:)); ylim([0 pi]);

% figure,
% plot(thetafinal(1,:),'.');







