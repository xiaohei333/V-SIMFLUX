function [mu,dmudtheta,PSF,PSFder] = poissonrate_new(params,theta,PupilMatrix,allzernikes,wavevector,wavevectorzimm)
% Returns the Poisson-rates for all pixels and all first order derivatives
% w.r.t. the parameters theta.

K = params.K;
% m = params.m;
Mx = params.Mx;
My = params.My;

fitmodel = params.fitmodel;
numparams = params.numparams;
% alpha = params.alpha;
% beta = params.beta;

switch fitmodel
        case 'xyz-azim-pola-diffusion'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        params.g2 = theta(8);
end

switch params.excitation
    case 'constant'
        P = 1;
        dPdazim = 0;
        dPdpola = 0;
    case 'zstack'
        P = ones(1,K);
        dPdazim = zeros(1,K);
        dPdpola = zeros(1,K);
end

% update pupil function
if contains(fitmodel,'aberrations')
    [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);
end

[FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
[PSF,PSFder] = get_psfs_derivatives(params,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);

% get Poisson rate and derivatives
mu = zeros(Mx,My,K);
dmudtheta = zeros(Mx,My,K,numparams);


for k = 1:K
    mu(:,:,k) = Nph*P(k)*PSF+Nbg/K;
         if strcmp(fitmodel,'xyz-azim-pola-diffusion')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = Nph*P(k)*PSFder(:,:,3);
            dmudtheta(:,:,k,4) = P(k)*PSF;
            dmudtheta(:,:,k,5) = 1/K;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,7) = Nph*P(k)*PSFder(:,:,5)+Nph*dPdpola(k)*PSF;
            dmudtheta(:,:,k,8) = Nph*P(k)*PSFder(:,:,6);
         end
end
end