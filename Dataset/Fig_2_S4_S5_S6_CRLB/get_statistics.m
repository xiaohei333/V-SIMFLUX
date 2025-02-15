function [err_std, crlb_mean] = get_statistics(object, est, crlb, outliers)
% calculation of standard deviation of the estiamtes

% outliers
outliers = outliers | (est(7,:)>deg2rad(160)) | (est(7,:)<deg2rad(20));

% error
err = est - object;

% convert orientations to half-sphere
flag = est(end-2,:)>pi;
est(end-2, flag) = est(end-2,flag)-pi;
est(end-1, flag) = pi-est(end-1,flag);

% calculate orientation error using unit vectors
Aest = Avec(est(end-2,:),est(end-1,:));
Aobj = Avec(object(end-2,:),object(end-1,:));

flip = sum(sign(Aobj)==sign(Aest))<2; % flipped match

norm_azim = mod(est(end-2,:)-object(end-2,:),pi);
absdiff_azim = min(pi-norm_azim,norm_azim);
absdiff_azim(norm_azim>pi/2) = -absdiff_azim(norm_azim>pi/2);

norm_pola = mod(est(end-1,:)-object(end-1,:),pi);
norm_pola(flip) = mod(pi-est(end-1,flip)-object(end-1,flip),pi);
absdiff_pola = min(pi-norm_pola,norm_pola);
absdiff_pola(norm_pola>pi/2) = -absdiff_pola(norm_pola>pi/2);

err(end-2,:) = absdiff_azim;
err(end-1,:) = absdiff_pola;

err_std = std(err(:, ~outliers), 0, 2);
crlb_mean = mean(crlb(:, ~outliers), 2);

end

% calculate spherical vector (A)
function [A] = Avec(azim,pola)

A = [cos(azim).*sin(pola); sin(azim).*sin(pola); cos(pola)];

end
