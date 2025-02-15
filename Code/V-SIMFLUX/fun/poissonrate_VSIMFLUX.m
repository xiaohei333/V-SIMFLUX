function [mu,dmudtheta] = poissonrate_VSIMFLUX(params,SIM,theta,PupilMatrix,allzernikes,wavevector,wavevectorzimm,ROIxy)

Mx = params.Mx;
My = params.My;
fitmodel = params.fitmodel;
numparams = params.numparams;

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

r = [theta(1)+ROIxy(2) theta(2)+ROIxy(1)];
g2 = theta(8);

L = SIM.L; % 条纹方向
K = SIM.K; % 相移步数
BetalStore = SIM.Betal; % a global angular offset
Chi = SIM.Chi; % size(Chi) = 1*L
lp = SIM.lp; % 条纹周期（nm）
m = SIM.m; % 调制深度
Phi = zeros(L,K);
P = zeros(L,K);


for l = 1:L
    for k = 1:K
        Betal = BetalStore(l);
        ql = (1/lp)*[cos(Betal) sin(Betal)];
        Psi = 2*pi*(k-1)/K + Chi(l);
        Phi(l,k) = 2*pi*ql*r' - Psi;
    end
end
for l = 1:L
    for k = 1:K
        P(l,k) = 1/(L*K)*( 1 + m(l)*cos(Phi(l,k)) );
    end
end

dipor = [sin(pola)*cos(azim) sin(pola)*sin(azim) cos(pola)];
AbsorbEff = zeros(1,L);
for l = 1:L
    Betal_E = BetalStore(l) + pi/2;
    AbsorbEff(l) = (1-g2)/2 + g2*([cos(Betal_E) sin(Betal_E) 0]*dipor')^2;
end

Eta = zeros(L,K);
dEtadtheta = zeros(L,K,5);
for l = 1:L
    for k = 1:K
        Eta(l,k) = P(l,k)*AbsorbEff(l) / ( (1-g2)/2 + sin(pola)^2*g2 /2 );
    end
end

Beta0 = (BetalStore(1) - pi/2 + BetalStore(2) - pi)/2;
%对x的导数
for l = 1:L
    for k = 1:K
        Betal = BetalStore(l);
        qlx = cos(Betal)/lp;
        dEtadtheta(l,k,1) = ( AbsorbEff(l)/((1-g2)/2+sin(pola)^2*g2/2) )*(-2)*pi*m(l)*qlx*sin(Phi(l,k))/(K*L);
    end
end
%对y的导数
for l = 1:L
    for k = 1:K
        Betal = BetalStore(l);
        qly = sin(Betal)/lp;
        dEtadtheta(l,k,2) = ( AbsorbEff(l)/((1-g2)/2+sin(pola)^2*g2/2) )*(-2)*pi*m(l)*qly*sin(Phi(l,k))/(K*L);
    end
end
%对dazim的导数
for l = 1:L
    for k = 1:K
        dEtadtheta(l,k,3) = P(l,k)*(-1)^l*( 4*g2*sin(pola)^2*sin(2*(Beta0 - azim)) )/( g2-2+g2*cos(2*pola) );
    end
end
%对dpola的导数
for l = 1:L
    for k = 1:K
        dEtadtheta(l,k,4) = P(l,k)*(-1)^l*( 4*(g2-1)*g2*sin(2*pola)*cos(2*(Beta0-azim)) )/( g2-2+g2*cos(2*pola) )^2;
    end
end
%对dg2的导数
for l = 1:L
    for k = 1:K
        dEtadtheta(l,k,5) = P(l,k)*(-1)^(l-1)*( 4*sin(pola)^2*cos(2*(Beta0-azim)) )/( g2-2+g2*cos(2*pola) )^2;
    end
end


% update pupil function
if contains(fitmodel,'aberrations')
    [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);
end

[FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
[PSF,PSFder] = get_psfs_derivatives(params,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);

% get Poisson rate and derivatives
mu = zeros(Mx,My,params.K,L,K);
dmudtheta = zeros(Mx,My,params.K,numparams,L,K);

for l = 1:L
    for k = 1:K
        mu(:,:,params.K,l,k) = Nph*Eta(l,k)*PSF+Nbg/(L*K);

        dmudtheta(:,:,params.K,1,l,k) = Nph.*( dEtadtheta(l,k,1).*PSF + Eta(l,k).*PSFder(:,:,1) );

        dmudtheta(:,:,params.K,2,l,k) = Nph.*( dEtadtheta(l,k,2).*PSF + Eta(l,k).*PSFder(:,:,2) );

        dmudtheta(:,:,params.K,3,l,k) = Nph.*( Eta(l,k).*PSFder(:,:,3) );

        dmudtheta(:,:,params.K,4,l,k) = Eta(l,k).*PSF;

        dmudtheta(:,:,params.K,5,l,k) = 1/(L*K);

        dmudtheta(:,:,params.K,6,l,k) = Nph.*( dEtadtheta(l,k,3).*PSF + Eta(l,k).*PSFder(:,:,4) );

        dmudtheta(:,:,params.K,7,l,k) = Nph.*( dEtadtheta(l,k,4).*PSF + Eta(l,k).*PSFder(:,:,5) );

        dmudtheta(:,:,params.K,8,l,k) = Nph.*( dEtadtheta(l,k,5).*PSF + Eta(l,k).*PSFder(:,:,6) );
    end
end

end