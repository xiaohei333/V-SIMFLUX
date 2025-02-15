close all;
clearvars;

% set parameters
params = set_parameters_vortex_sim;
params.Ncfg = 1000; % generate and fit Ncfg randomized model PSFs.
params.flg_parallel = 1;

[wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);

object = zeros(params.numparams,params.Ncfg);
allpsfs = zeros(params.Mx,params.My,params.K,params.Ncfg);
dallpsfsdtheta = zeros(params.Mx,params.My,params.K,params.numparams,params.Ncfg);

g2store = [0.001,0.25:0.25:0.75,0.999];
azimstore =  deg2rad(0:10:180);
polastore =  deg2rad(20:10:90);

CrlbMeanStore_Vortex = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2));
ThetaMeanStore_Vortex = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2));
ThetaStdStore_Vortex = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2));
OutliersMeanStore_Vortex = zeros(size(azimstore,2),size(polastore,2),size(g2store,2));


ObjectStore_Vortex = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);
ThetaStore_Vortex = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);
CrlbStore_Vortex = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);
OutliersStore_Vortex = zeros(size(azimstore,2),size(polastore,2),size(g2store,2),params.Ncfg);

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

                object(:,jj) = [dx dy dz Nphotons Nbackground dazim dpola dg2];

                [allpsfs(:,:,params.K,jj) , dallpsfsdtheta(:,:,params.K,:,jj)] ...
                    = poissonrate(params,object(:,jj),PupilMatrix,allzernikes,wavevector,wavevectorzimm);

            end
            mu_truth = allpsfs;
            dmudtheta_truth = dallpsfsdtheta;
            allspots = 1e12*imnoise(allpsfs*1e-12,'poisson');

            % 进行拟合
            [thetainit] = initialvalues(allspots,params);
            [thetastore,mu,dmudtheta,merit,numiters] = localization(allspots,thetainit,params);
            theta = thetastore(:,:,end);

            [crlb,rcondstore] = get_fisher_crlb(params,mu_truth,dmudtheta_truth);
            [outliers] = get_outliers(theta,merit,numiters,params);
            [thetafinal,thetamean,thetastd,crlbmean] = get_statistics(params,object,theta,crlb,outliers);

            ThetaMeanStore_Vortex(:,mm,nn,ii) = thetamean;
            ThetaStdStore_Vortex(:,mm,nn,ii) = thetastd;
            CrlbMeanStore_Vortex(:,mm,nn,ii) = crlbmean;
            OutliersMeanStore_Vortex(mm,nn,ii) = sum(outliers,'all');

            ObjectStore_Vortex(:,mm,nn,ii,:) = object;
            ThetaStore_Vortex(:,mm,nn,ii,:) = theta;
            CrlbStore_Vortex(:,mm,nn,ii,:) = crlb;
            OutliersStore_Vortex(mm,nn,ii,:) = outliers;

        end
        fprintf('目前进程为 %.2f %% \n',100*ii/size(g2store,2));
        fprintf('子进程为 %.2f %% \n',100*((mm-1)*size(polastore,2)+nn)/size(azimstore,2)/size(polastore,2));
    end
    % fprintf('Process is %.2f %% \n',100*ii/size(SBRstore,2));
end

%% 结果展示

% save data.mat;
save('Vortex.mat','ThetaMeanStore_Vortex','ThetaStdStore_Vortex','CrlbMeanStore_Vortex','OutliersMeanStore_Vortex',...
    'ObjectStore_Vortex','ThetaStore_Vortex','CrlbStore_Vortex','OutliersStore_Vortex');


figure,
for ii = 1:params.numparams
subplot(2,4,ii);
plot(g2store,squeeze(mean(CrlbMeanStore_Vortex(ii,:,:,:),[2 3]))); hold on;
plot(g2store,squeeze(mean(ThetaStdStore_Vortex(ii,:,:,:),[2 3]))); 
legend('CRLB','Precision');
end

figure,
for ii = 1:params.numparams
subplot(2,4,ii);
plot(g2store,squeeze(mean(ThetaMeanStore_Vortex(ii,:,:,:),[2 3])));
end

figure,
for ii = 1:params.numparams
subplot(2,4,ii);
plot(g2store,squeeze(mean(OutliersMeanStore_Vortex,[1 2])),'*');
end

A = squeeze(CrlbStore_Vortex(6,:,:,4));
B = rad2deg(mean(A(:,21:91),2));
C = rad2deg(mean(A(1:10:181,21:10:91),2));

figure,
imagesc(A); axis image;
figure,
plot(1:181,B); hold on;
plot(1:10:181,C);
ylim([0 15]);
