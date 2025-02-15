clear;
load("spl_images.mat");

n = 1;
dphi = 0.001;

%% DISPLAY
figure;
for jj = 1:3
    for mm = 1:size(images, 5)
        % normalization for display
        dispMin = min(images(:,:,:,:,mm), [], 'all');
        dispMax = max(images(:,:,:,:,mm), [], 'all');

        for nn = 1:size(images, 3)
            subplot(size(images, 5), 6, 6*(mm-1)+3*(nn-1)+jj);
            imagesc(images(:,:,nn,jj,mm));
            axis image off; colormap(gca, 'hot');
            clim([dispMin dispMax]);
        end
    end
end
