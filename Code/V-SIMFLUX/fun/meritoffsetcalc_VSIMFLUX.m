function meritoffset = meritoffsetcalc_VSIMFLUX(allspots,varfit)
% This function calculates the merit function offset for an accurate
% determination of the log-likelihood.
%

[Mx,My,PK,L,K,Ncfg] = size(allspots);
meritoffset = zeros(Ncfg,1,'Like',allspots);
meritoffsetTemp = zeros(Ncfg,1,L,K,'Like',allspots);

for jcfg = 1:Ncfg
 for l=1:L
  for k=1:K
    dummat = allspots(:,:,:,l,k,jcfg);
    % set negative and zero pixels values to one/10 to avoid log-singularity
    dummat = max(dummat,ones(size(dummat))/10);
    meritoffsetTemp(jcfg,1,l,k) = 0;
    for ii=1:Mx
        for jj = 1:My
            for kk = 1:PK
                meritoffsetTemp(jcfg,1,l,k) = meritoffsetTemp(jcfg,1,l,k)-gammln(dummat(ii,jj,kk)+1+varfit);
            end
        end
    end
  end
 end
end

meritoffset = sum(meritoffsetTemp,[3 4]);

end
