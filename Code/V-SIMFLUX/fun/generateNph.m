function Nbackground = generateNph(SIM,theta,ROIxy)

Nph = theta(4);
azim = theta(6);
pola = theta(7);
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
for l = 1:L
    for k = 1:K
        Eta(l,k) = P(l,k)*AbsorbEff(l) / ( (1-g2)/2 + sin(pola)^2*g2 /2 );
    end
end

Nbackground = Nph.*Eta;
Nbackground = Nbackground';
Nbackground = Nbackground(:);

end