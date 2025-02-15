close all;
clearvars;

% set parameters
params = set_parameters_vortex_sim;
params.Ncfg = 1; % generate and fit Ncfg randomized model PSFs.
params.flg_parallel = 1;

[wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);

object = zeros(params.numparams,params.Ncfg);
allpsfs = zeros(params.Mx,params.My,params.K,params.Ncfg);
dallpsfsdtheta = zeros(params.Mx,params.My,params.K,params.numparams,params.Ncfg);

g2store = [0.001,0.25:0.25:0.75,0.999];
azimstore =  deg2rad(0:180);
polastore =  deg2rad(0:180);

CrlbStore_Vortex = zeros(params.numparams,size(azimstore,2),size(polastore,2),size(g2store,2));

for ii = 1:size(g2store,2)
    for mm = 1:size(azimstore,2)
        for nn = 1:size(polastore,2)

            % 产生图片
            for jj = 1:params.Ncfg

                % true parameters
                dx = 0;
                dy = 0;
                dz = 0;
                Nphotons = 4000;
                Nbackground = 10;
                dazim = azimstore(mm);
                dpola = polastore(nn);
                dg2 = g2store(ii);

                object(:,jj) = [dx dy dz Nphotons Nbackground dazim dpola dg2];
                [allpsfs(:,:,params.K,jj) , dallpsfsdtheta(:,:,:,:,jj)] ...
                    = poissonrate(params,object(:,jj),PupilMatrix,allzernikes,wavevector,wavevectorzimm);


            end

            [CRLBstore,rcondstore] = get_fisher_crlb(params,allpsfs,dallpsfsdtheta);
            CrlbStore_Vortex(:,mm,nn,ii) = mean(CRLBstore,2);


        end
        fprintf('目前进程为 %.2f %% \n',100*ii/size(g2store,2));
        fprintf('子进程为 %.2f %% \n',100*((mm-1)*size(polastore,2)+nn)/size(azimstore,2)/size(polastore,2));
    end
    % fprintf('Process is %.2f %% \n',100*ii/size(SBRstore,2));
end

%% 展示结果

figure,
imagesc(rad2deg(squeeze(CrlbStore_Vortex(6,:,:,5))));
clim([0 15]);

figure,
imagesc(squeeze(CrlbStore_Vortex(1,:,:,5)));