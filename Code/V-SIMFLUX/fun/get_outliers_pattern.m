function [m1,m2,outliers_pattern1,outliers_pattern2] = get_outliers_pattern(amp1, offset1, amp2, offset2, SIM_m, N, bg)

L = 2;
K = 3;
SBR_threshold = 100;
flag1 = sum( N(:,1:3) <= SBR_threshold*bg/L/K , 2) < 2;
flag2 = sum( N(:,4:6) <= SBR_threshold*bg/L/K , 2) < 2;
m1 = amp1(flag1)./offset1(flag1);
m2 = amp2(flag2)./offset2(flag2);

% m1 = amp1./offset1;
% m2 = amp2./offset2;
m1_mean = mean(m1);
m1_std = std(m1);
m2_mean = mean(m2);
m2_std = std(m2);

outliers_pattern1 = (abs(m1)>SIM_m(1)+0.15)|(abs(m1)<SIM_m(1)-0.15);
outliers_pattern1 = outliers_pattern1|(m1<m1_mean-3*m1_std)|(m1>m1_mean+3*m1_std);

outliers_pattern2 = (abs(m2)>SIM_m(2)+0.15)|(abs(m2)<SIM_m(2)-0.15);
outliers_pattern2 = outliers_pattern2|(m2<m2_mean-3*m2_std)|(m2>m2_mean+3*m2_std);

end