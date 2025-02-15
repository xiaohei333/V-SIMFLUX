function [ m, outliers_pattern ] = modulation_fitting(xcoor, ycoor, N, SIM_est)

K = SIM_est.K;
Betal = SIM_est.Betal;
Chi = SIM_est.Chi;
lp = SIM_est.lp;
r = [xcoor ycoor];

N_1 =  N(:,1:3);
flag = sum( N_1 <= 100*10/6 , 2) < 2;
P_1 = N_1(flag,:)./sum(N_1(flag,:),2);


m(1) = lsqcurvefit(@(m,r) IntensityProfile(m,r,K,Betal(1),Chi(1),lp),0.9,r(flag,:),P_1);

end


function P = IntensityProfile(m,r,K,Betal,Chi,lp)

Phi = zeros(size(r,1),K);
ql = (1/lp)*[cos(Betal) sin(Betal)];
for k = 1:K
    Psi = 2*pi*(k-1)/K + Chi;
    Phi(:,k) = 2*pi*ql*r' - Psi;
end

P = 1 + m*cos(Phi);
P = P./sum(P,2);

end