clear;
load("images.mat");

n = 1;
dphi = 0.001;
blackboarder = 20;

%% DISPLAY
figure;
for jj = 1:3
    % generating cosine fringes
    phi = (-1:dphi:1)*n*2*pi - 2*(jj-1)*pi/3;
    i = cos(phi./2).^2;

    I = repmat(i, [2001 1]);
    It = repmat(i',[1 2001]);

    % insert black border
    I1 = zeros(2*blackboarder+2001, 2*blackboarder+2001);
    I2 = I1;
    I1(blackboarder+1:blackboarder+2001, blackboarder+1:blackboarder+2001) = I;
    I2(blackboarder+1:blackboarder+2001, blackboarder+1:blackboarder+2001) = It;
    
    % display fringes
    for kk = 1:2
        if kk == 1
            subplot(1+length(images), 6, jj);
            imagesc(I1);
        else
            subplot(1+length(images), 6, jj+3);
            imagesc(I2);
        end
        axis image off; colormap(gca, 'gray');
        hold on;
        plot(1001, 1001, 'p', 'LineWidth', 1,...
            'MarkerEdgeColor', [245 128 42]./255,...
            'MarkerFaceColor', [245 128 42]./255,...
            'MarkerSize', 5);
        hold off;
    end
    
    for mm = 1:length(images)
        % normalization for display
        dispMin = min(images(mm).psf, [], 'all');
        dispMax = max(images(mm).psf, [], 'all');

        for nn = 1:size(images(mm).psf, 3)
            subplot(1+length(images), 6, 6*mm+3*(nn-1)+jj);
            imagesc(images(mm).psf(:,:,nn,jj));
            axis image off; colormap(gca, 'hot');
            clim([dispMin dispMax]);
        end
    end

end

