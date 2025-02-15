clear;
load("psf.mat");

%% DISPLAY
figure;

for i = 1:size(psf, 4)
    % normalization for display
    dispMin = min(psf(:,:,:,i), [], 'all');
    dispMax = max(psf(:,:,:,i), [], 'all');

    for j = 1:5
        subplot(3, 5, (i-1)*5+j);
        imagesc(psf(:,:,j,i));
        axis image off; colormap(gca, 'hot');
        clim([dispMin dispMax]);
    end

end
