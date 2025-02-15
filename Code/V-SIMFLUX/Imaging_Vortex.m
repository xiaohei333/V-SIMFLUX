close all;
clearvars;

localNum = 1;
azimstore = pi/5*(1:10);
polastore = pi/2;
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

params = set_parameters_vortex_sim;
params.Ncfg = localNum*size(azimstore,2)*size(polastore,2); % generate and fit Ncfg randomized model PSFs.
params.flg_parallel = 1;


% figure;
% plot3(xstore(1:10),ystore(1:10),zstore(1:10),'*');
% xlabel('X axis(nm)');
% ylabel('Y axis(nm)');
% zlabel('Z axis(nm)');
% xlim([-200 200]);
% ylim([-200 200]);

%% generate ground truth PSFs
[wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);

allpsfs = zeros(params.Mx,params.My,params.K,params.Ncfg);
object = zeros(params.numparams,params.Ncfg);

for ii = 1:size(azimstore,2)
    for jj = 1:size(polastore,2)
        for mm = 1:localNum

            index = (ii-1)*size(polastore,2)*localNum + (jj-1)*localNum + mm;
            % true parameters
            dx = xstore(ii,jj);
            dy = ystore(ii,jj);
            dz = zstore(ii,jj);
            Nphotons = 4000;
            Nbackground = 10;
            dazim = azimstore(ii);
            dpola = polastore(jj);
            dg2 = 0;

            % generate object and PSFs
            object(:,index) = [dx dy dz Nphotons Nbackground dazim dpola dg2];
            [allpsfs(:,:,params.K,index) dallpsfsdtheta(:,:,params.K,:,index)] ...
                = poissonrate(params,object(:,index),PupilMatrix,allzernikes,wavevector,wavevectorzimm);

        end
    end
end

% add noise
allspots = 1e12*imnoise(allpsfs*1e-12,'poisson');

%% MLE fitting routine
[thetainit] = initialvalues(allspots,params);
[thetastore,mu,dmudtheta,merit,numiters] = localization(allspots,thetainit,params);
theta = thetastore(:,:,end);


[crlb,rcondstore] = get_fisher_crlb(params,mu,dmudtheta);
[outliers] = get_outliers(theta,merit,numiters,params);
[thetafinal,thetamean,thetastd,crlbmean] = get_statistics(params,object,theta,crlb,outliers);

%% 效果展示

x_coor = theta(1,:);
y_coor = theta(2,:);
figure,
plot(x_coor(~outliers),y_coor(~outliers),'.');
xlim([-80 80]);
ylim([-80 80]);

% tmpcfg = theta(end-2,:)>pi;
% theta(end-2,tmpcfg) = theta(end-2,tmpcfg)-pi;
% theta(end-1,tmpcfg) = pi-theta(end-1,tmpcfg);
figure,
subplot(1,2,1);plot(theta(6,~outliers));
subplot(1,2,2);plot(theta(7,~outliers));


figure,
subplot(1,3,1);
histogram(thetafinal(1,~outliers)); title('X与真值的差（去除异常点）');
subplot(1,3,2);
histogram(thetafinal(2,~outliers)); title('Y与真值的差（去除异常点）');
subplot(1,3,3);
histogram(thetafinal(6,~outliers)); title('azim与真值的差（去除异常点）');










