function [lp, resnormx] = FitPhasePlane_SIMFLUX(x,y,p2, lp0, fast)

 if nargin<5
        fast = 0;
 end
    
 if fast==0
     options = optimset('Display','off','MaxFunEvals',1000,'MaxIter',100,'TolFun',1e-5,'LargeScale','on');
 else
     options = optimset('Display','off','MaxFunEvals',100,'MaxIter',5,'TolFun',1e-4,'LargeScale','off'); %fast==1时候
 end

[lpx,resnormx,~,~]=lsqcurvefit(@(xp, xdata)CalPhaseCos(xp, x, y), ...
                   [lp0(1) lp0(2) 0],x,[sin(p2) cos(p2)],[],[],options);

lp = [lpx(1), lpx(2), checkPhase(lpx(3))];

end