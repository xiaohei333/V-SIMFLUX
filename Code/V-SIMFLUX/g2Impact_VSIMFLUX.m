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

% g2store = 0:0.25:1;
g2store = [0.001,0.25:0.25:0.75,0.999];
azimstore =  deg2rad(0:10:180);
polastore =  deg2rad(20:10:90);

% CrlbMeanStore_V = zeros(params.numparams,size(SBRstore,2));
% ThetaMeanStore_V = zeros(params.numparams,size(SBRstore,2));
% ThetaStdStore_V = zeros(params.numparams,size(SBRstore,2));
CrlbMeanStore_VSIMFLUX = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2));
ThetaMeanStore_VSIMFLUX = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2));
ThetaStdStore_VSIMFLUX = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2));
OutliersMeanStore_VSIMFLUX = zeros(size(azimstore,2),size(polastore,2),size(g2store,2));


ObjectStore_VSIMFLUX = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);
ThetaStore_VSIMFLUX = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);
CrlbStore_VSIMFLUX = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);
OutliersStore_VSIMFLUX = zeros(size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);


%% 开始计算

for ii = 1:size(g2store,2)
    for mm = 1:size(azimstore,2)
        for nn = 1:size(polastore,2)

            % 产生图片
            for jj = 1:params.Ncfg

                % true parameters
                dx = (1-2*rand)*params.pixelsize;
                dy = (1-2*rand)*params.pixelsize;
                dz = 0;
                Nphotons = 4000;
                Nbackground = 10;
                dazim = azimstore(mm);
                dpola = polastore(nn);
                dg2 = g2store(ii);

                ROIxy = [0 0];

                object(:,jj) = [dx dy dz Nphotons Nbackground dazim dpola dg2];

                [allpsfs(:,:,params.K,:,:,jj) , dallpsfsdtheta(:,:,:,:,:,:,jj)] ...
                    = poissonrate_VSIMFLUX(params,SIM,object(:,jj),PupilMatrix,allzernikes,wavevector,wavevectorzimm,ROIxy);

            end
            mu_truth = allpsfs;
            dmudtheta_truth = dallpsfsdtheta;
            allpsfs = 1e12*imnoise(allpsfs*1e-12,'poisson');

            % 进行拟合
            allpsfs_V(:,:,1,:) = squeeze(sum(allpsfs,[4 5]));
            [thetainit_V] = initialvalues(allpsfs_V,params);
            [thetastore_V,mu,dmudtheta,merit,numiters] = localization(allpsfs_V,thetainit_V,params);
            theta_V = thetastore_V(:,:,end);
            thetainit_VSIMFLUX = theta_V;

            roixy0 = zeros(params.Ncfg,2);
            [thetastore_VSIMFLUX,mu_VSIMFLUX,dmudtheta_VSIMFLUX,merit_VSIMFLUX,numiters_VSIMFLUX] ...
                = localization_VSIMFLUX(allpsfs,thetainit_VSIMFLUX,params,SIM,roixy0);
            theta_VSIMFLUX = thetastore_VSIMFLUX(:,:,end);

            % 数据分析
            % [crlb,rcondstore] = get_fisher_crlb(params,mu,dmudtheta);
            % [outliers] = get_outliers(theta_V,merit,numiters,params);
            % [thetafinal,thetamean,thetastd,crlbmean] = get_statistics(params,object,theta_V,crlb,outliers);

            % [crlb_VSIMFLUX,rcondstore_VSIMFLUX] = get_fisher_crlb_VSIMFLUX(params,SIM,mu_VSIMFLUX,dmudtheta_VSIMFLUX);
            [crlb_VSIMFLUX,rcondstore_VSIMFLUX] = get_fisher_crlb_VSIMFLUX(params,SIM,mu_truth,dmudtheta_truth);
            [outliers_VSIMFLUX] = get_outliers(theta_VSIMFLUX,merit_VSIMFLUX,numiters_VSIMFLUX,params);
            [thetafinal_VSIMFLUX,thetamean_VSIMFLUX,thetastd_VSIMFLUX,crlbmean_VSIMFLUX] ...
                = get_statistics(params,object,theta_VSIMFLUX,crlb_VSIMFLUX,outliers_VSIMFLUX);


            % ThetaMeanStore_V(:,ii) = thetamean;
            % ThetaStdStore_V(:,ii) = thetastd;
            % CrlbMeanStore_V(:,ii) = crlbmean;

            ThetaMeanStore_VSIMFLUX(:,mm,nn,ii) = thetamean_VSIMFLUX;
            ThetaStdStore_VSIMFLUX(:,mm,nn,ii) = thetastd_VSIMFLUX;
            CrlbMeanStore_VSIMFLUX(:,mm,nn,ii) = crlbmean_VSIMFLUX;
            OutliersMeanStore_VSIMFLUX(mm,nn,ii) = sum(outliers_VSIMFLUX,'all');

            ObjectStore_VSIMFLUX(:,mm,nn,ii,:) = object;
            ThetaStore_VSIMFLUX(:,mm,nn,ii,:) = theta_VSIMFLUX;
            CrlbStore_VSIMFLUX(:,mm,nn,ii,:) = crlb_VSIMFLUX;
            OutliersStore_VSIMFLUX(mm,nn,ii,:) = outliers_VSIMFLUX;

        end
        fprintf('目前进程为 %.2f %% \n',100*ii/size(g2store,2));
        fprintf('子进程为 %.2f %% \n',100*((mm-1)*size(polastore,2)+nn)/size(azimstore,2)/size(polastore,2));
    end
    % fprintf('Process is %.2f %% \n',100*ii/size(SBRstore,2));
end

%% 结果展示

% save data.mat;

save('VSIMFLUX.mat','ThetaMeanStore_VSIMFLUX','ThetaStdStore_VSIMFLUX','CrlbMeanStore_VSIMFLUX','OutliersMeanStore_VSIMFLUX',...
    'ObjectStore_VSIMFLUX','ThetaStore_VSIMFLUX','CrlbStore_VSIMFLUX','OutliersStore_VSIMFLUX');

figure,
for ii = 1:params.numparams
subplot(2,4,ii);
plot(g2store,squeeze(mean(CrlbMeanStore_VSIMFLUX(ii,:,:,:),[2 3]))); hold on;
plot(g2store,squeeze(mean(ThetaStdStore_VSIMFLUX(ii,:,:,:),[2 3]))); 
legend('CRLB','Precision');
end

figure,
for ii = 1:params.numparams
subplot(2,4,ii);
plot(g2store,squeeze(mean(ThetaMeanStore_VSIMFLUX(ii,:,:,:),[2 3])));
end

figure,
for ii = 1:params.numparams
subplot(2,4,ii);
plot(g2store,squeeze(mean(OutliersMeanStore_VSIMFLUX,[1 2])),'*');
end

% figure,
% subplot(1,2,1);
% imagesc(squeeze(CrlbMeanStore_VSIMFLUX(6,:,:,5)));
% subplot(1,2,2);
% imagesc(squeeze(ThetaStdStore_VSIMFLUX(6,:,:,5)));


figure,
subplot(1,2,1);
imagesc(squeeze(CrlbMeanStore_VSIMFLUX(3,:,:,2)));
subplot(1,2,2);
imagesc(squeeze(ThetaStdStore_VSIMFLUX(3,:,:,2)));
