close all;
clearvars;

localNum = 100;
azimstore = pi/5*(1:10);
polastore = pi/2;
% localNum = 1;
% azimstore = pi/4;
% polastore = pi/4;
R = 50;
xstore = zeros(size(azimstore,2),size(polastore,2));
ystore = zeros(size(xstore));
zstore = zeros(size(xstore));

for ii = 1:size(azimstore,2)
    for jj = 1:size(polastore,2) 
        xstore(ii,jj) = R*sin(polastore(jj))*cos(azimstore(ii));
        ystore(ii,jj) = R*sin(polastore(jj))*sin(azimstore(ii));
        zstore(ii,jj) = R*cos(polastore(jj));
    end
end

% figure;
% plot3(xstore(1:10),ystore(1:10),zstore(1:10),'*');
% xlabel('X axis(nm)');
% ylabel('Y axis(nm)');
% zlabel('Z axis(nm)');
% xlim([-200 200]);
% ylim([-200 200]);

params = set_parameters_vortex_sim;
params.Ncfg = localNum*size(azimstore,2)*size(polastore,2); % generate and fit Ncfg randomized model PSFs.
params.flg_parallel = 1;

SIM.L = 2; % 条纹方向
SIM.K = 3; % 相移步数
SIM.Betal = [pi/2 pi]; % a global angular offset
SIM.Chi = [0 0]; % size(Chi) = 1*L
SIM.lp = 561/2/1.49; % 条纹周期（nm）
SIM.m = [0.9 0.9]; % 调制深度

[wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);

object = zeros(params.numparams,params.Ncfg);
allpsfs = zeros(params.Mx,params.My,params.K,SIM.L,SIM.K,params.Ncfg);
dallpsfsdtheta = zeros(params.Mx,params.My,params.K,params.numparams,SIM.L,SIM.K,params.Ncfg);

%% generate ground truth PSFs
for ii = 1:size(azimstore,2)
    for jj = 1:size(polastore,2)
        for mm = 1:localNum

            index = (ii-1)*size(polastore,2)*localNum + (jj-1)*localNum + mm;
            % true parameters
            dx = xstore(ii,jj);
            dy = ystore(ii,jj);
            dz = zstore(ii,jj);
            Nphotons = 2000;
            Nbackground = 10;
            dazim = azimstore(ii);
            dpola = polastore(jj);
            dg2 = 0.75;

            ROIxy = [0 0];
            object(:,index) = [dx dy dz Nphotons Nbackground dazim dpola dg2];

            [allpsfs(:,:,params.K,:,:,index) , dallpsfsdtheta(:,:,:,:,:,:,index)] ...
                = poissonrate_VSIMFLUX(params,SIM,object(:,index),PupilMatrix,allzernikes,wavevector,wavevectorzimm,ROIxy);

        end
    end
end

% add noise
allspots = 1e12*imnoise(allpsfs*1e-12,'poisson');

% figure,
% imagesc(squeeze(sum(allpsfs,[4 5])));
% 
figure,
for l = 1:2
    for k = 1:3
        subplot(2,3,(l-1)*3+k);
        imagesc(squeeze(allspots(:,:,1,l,k))); clim([0 15]); axis image;
    end
end


%% fitting parameters of illumination pattern

% 使用求和后的总图，初步估计各个参数
allspots_V(:,:,params.K,:) = squeeze( sum(allspots,[4 5]) );
[thetainit_V] = initialvalues(allspots_V,params);
[thetastore_V,mu,dmudtheta,merit,numiters] = localization(allspots_V,thetainit_V,params);
theta_V = thetastore_V(:,:,end);
% figure,
% plot(theta_V(1,:),theta_V(2,:),'.');
% xlim([-80 80]);
% ylim([-80 80]);

% 将VSIMFLUX的图片排列成一列
allspots_VSIMFLUX = zeros(params.Mx,params.My,params.K,SIM.L*SIM.K*params.Ncfg);
for ii = 1:params.Ncfg
    for l = 1:SIM.L
        for k = 1:SIM.K
            index = SIM.L*SIM.K*(ii-1)+SIM.K*(l-1)+k;
            allspots_VSIMFLUX(:,:,params.K,index) = allspots(:,:,params.K,l,k,ii);
        end
    end
end
% figure,
% subplot(1,3,1);imagesc(squeeze(allspots_VSIMFLUX(:,:,1,1))); clim([0 15]);
% subplot(1,3,2);imagesc(squeeze(allspots_VSIMFLUX(:,:,1,2))); clim([0 15]);
% subplot(1,3,3);imagesc(squeeze(allspots_VSIMFLUX(:,:,1,3))); clim([0 15]);
% figure,
% subplot(1,3,1);imagesc(squeeze(allspots(:,:,1,1,1,1))); clim([0 15]);
% subplot(1,3,2);imagesc(squeeze(allspots(:,:,1,1,2,1))); clim([0 15]);
% subplot(1,3,3);imagesc(squeeze(allspots(:,:,1,1,3,1))); clim([0 15]);


%%
% 估计各个子图的光子数
params.Ncfg = size(allspots_VSIMFLUX,4);
thetainit_VSIMFLUX = zeros(params.numparams,params.Ncfg);
% 将之前总图估计得到的 x,y,z,Nb/6,phi,theta,g2 替换thetainitSIMFLUX中的初始值
for ii = 1:(params.Ncfg/SIM.L/SIM.K)
    index = (ii-1)*SIM.L*SIM.K + (1:SIM.L*SIM.K);
    thetainit_VSIMFLUX(1:3,index) = repmat(theta_V(1:3,ii),[1 6]);
    thetainit_VSIMFLUX(5,index) = repmat(theta_V(5,ii),[1 6])./(SIM.L*SIM.K);
    thetainit_VSIMFLUX(6:8,index) = repmat(theta_V(6:8,ii),[1 6]);
end
% 估计Nphotons的初始值
for ii = 1:params.Ncfg
    bg = thetainit_VSIMFLUX(5,ii);
    thetainit_VSIMFLUX(4,ii) = 2.5*sum(max(allspots_VSIMFLUX(:,:,:,ii)-bg,1e-12),'all'); % photonflux = 2.5
end

% 只对Nphotons这一参数进行最优化
[Nstore_VSIMFLUX,muSIMFLUX,dmudthetaSIMFLUX,meritSIMFLUX,numitersSIMFLUX] ...
    = localization_N(allspots_VSIMFLUX,thetainit_VSIMFLUX,params);
N_VSIMFLUX = Nstore_VSIMFLUX(:,:,end);

N = zeros(params.Ncfg/SIM.L/SIM.K,SIM.L*SIM.K);
N(:,1) = N_VSIMFLUX(1:6:end);
N(:,2) = N_VSIMFLUX(2:6:end);
N(:,3) = N_VSIMFLUX(3:6:end);
N(:,4) = N_VSIMFLUX(4:6:end);
N(:,5) = N_VSIMFLUX(5:6:end);
N(:,6) = N_VSIMFLUX(6:6:end);

%%
%估计条纹的ml和相位
phase0 = 2*pi/3;
[phase1, amp1, offset1] = CalPhase(N(:,1:3), phase0);
[phase2, amp2, offset2] = CalPhase(N(:,4:6), phase0);
[m1,m2,outliers_pattern1,outliers_pattern2] = get_outliers_pattern(amp1, offset1, amp2, offset2, SIM.m, N, theta_V(5,:)');

%估计条纹的方向，周期，初始相位
lp1 = SearchPhasePlane(theta_V(1,:)',theta_V(2,:)',phase1,[85/180*pi 95/180*pi]);
lp2 = SearchPhasePlane(theta_V(1,:)',theta_V(2,:)',phase2,[175/180*pi 185/180*pi]);

SIM_est.L = SIM.L; % 条纹方向
SIM_est.K = SIM.K; % 相移步数
SIM_est.Betal = [lp1(1) lp2(1)]; % a global angular offset
SIM_est.Chi = [lp1(3) lp2(3)]; % size(Chi) = 1*L
SIM_est.lp = ( 2*pi/lp1(2) + 2*pi/lp2(2) )/2; % 条纹周期（nm）
SIM_est.m = [mean(m1(~outliers_pattern1)) mean(m2(~outliers_pattern2))];



% figure,
% subplot(1,2,1),plot(m1);
% subplot(1,2,2),plot(m2);
% figure,
% subplot(1,2,1),plot(m1(~outliers_pattern1));
% subplot(1,2,2),plot(m2(~outliers_pattern2));

%% 开始VSIMFLUX拟合

params.Ncfg = localNum*size(azimstore,2)*size(polastore,2);

thetainit_VSIMFLUX = theta_V;
roixy0 = zeros(params.Ncfg,2);
[thetastore_VSIMFLUX,mu_VSIMFLUX,dmudtheta_VSIMFLUX,merit_VSIMFLUX,numiters_VSIMFLUX] ...
    = localization_VSIMFLUX(allspots,thetainit_VSIMFLUX,params,SIM_est,roixy0);
theta_VSIMFLUX = thetastore_VSIMFLUX(:,:,end);

[crlb_VSIMFLUX,rcondstore_VSIMFLUX] = get_fisher_crlb_VSIMFLUX(params,SIM,allpsfs,dallpsfsdtheta);
[outliers_VSIMFLUX] = get_outliers(theta_VSIMFLUX,merit_VSIMFLUX,numiters_VSIMFLUX,params);
[thetafinal_VSIMFLUX,thetamean_VSIMFLUX,thetastd_VSIMFLUX,crlbmean_VSIMFLUX] ...
    = get_statistics(params,object,theta_VSIMFLUX,crlb_VSIMFLUX,outliers_VSIMFLUX);

%% 效果展示

x_coor = theta_VSIMFLUX(1,:);
y_coor = theta_VSIMFLUX(2,:);
figure,
plot(x_coor(~outliers_VSIMFLUX),y_coor(~outliers_VSIMFLUX),'.');
xlim([-80 80]);
ylim([-80 80]);

% tmpcfg = theta_VSIMFLUX(end-2,:)>pi;
% theta_VSIMFLUX(end-2,tmpcfg) = theta_VSIMFLUX(end-2,tmpcfg)-pi;
% theta_VSIMFLUX(end-1,tmpcfg) = pi-theta_VSIMFLUX(end-1,tmpcfg);
figure,
subplot(1,2,1);plot(theta_VSIMFLUX(6,~outliers_VSIMFLUX));
subplot(1,2,2);plot(theta_VSIMFLUX(7,~outliers_VSIMFLUX));


figure,
subplot(1,3,1);
histogram(thetafinal_VSIMFLUX(1,~outliers_VSIMFLUX)); title('X与真值的差（去除异常点）');
subplot(1,3,2);
histogram(thetafinal_VSIMFLUX(2,~outliers_VSIMFLUX)); title('Y与真值的差（去除异常点）');
subplot(1,3,3);
histogram(thetafinal_VSIMFLUX(6,~outliers_VSIMFLUX)); title('azim与真值的差（去除异常点）');
