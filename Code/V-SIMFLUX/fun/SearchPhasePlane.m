function [lp_optimized] = SearchPhasePlane(x,y,p2,searchrange)

    xs = x; 
    ys = y;
    p2s = p2;
    searchstart = min(searchrange);
    searchend = max(searchrange);  

    lplist1 = searchstart:0.005:searchend;
    lplist2 = 2*pi./(230:-1:210); %条纹周期从210nm搜索到230nm，步长为1nm

    [aa, bb] = meshgrid(lplist1, lplist2);
    resnormmap = zeros(size(aa)); 
    lplist = zeros(length(aa(:)), 3);

    for m=1:length(aa(:))
       [lp, resnormx] = FitPhasePlane_SIMFLUX(xs,ys,p2s, [aa(m) bb(m) 0], 0);
       resnormmap(m) = resnormx;  % 即在残差的平方和 16*21double
       lplist(m,:) = lp; 
    end

    mask = find(resnormmap <(min(resnormmap(:))+1));
    lplist = lplist(mask, :);
    options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-7,'LargeScale','on');
    [lpx,resnormx,~,~]=lsqcurvefit(@(xp, xdata)CalPhaseCos(xp, x, y), ...
    [mean(lplist(:,1)) mean(lplist(:,2)) 0],x,[sin(p2) cos(p2)],[],[],options);
    lp_optimized = [lpx(1), lpx(2), checkPhase(lpx(3))];

end