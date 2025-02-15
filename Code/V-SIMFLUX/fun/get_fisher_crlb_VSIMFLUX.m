function [CRLBstore,rcondstore] = get_fisher_crlb_VSIMFLUX(params,SIM,mustore,dmudthetastore)
% fprintf('\nFisher information: '); tic;
% This function calculates the Fisher-matrix and the Cramer-Rao Lower Bound
% for the parameters found.

keps = 1e3*eps;
Ncfg = params.Ncfg;
numparams = params.numparams;

L = SIM.L;
K = SIM.K;

CRLBstore = zeros(numparams,Ncfg);
rcondstore = zeros(1,Ncfg);
for jcfg = 1:Ncfg

    FisherTemp = zeros(numparams,numparams,L,K);
    for l = 1:L
        for k = 1:K

            mu = mustore(:,:,params.K,l,k,jcfg);
            dmudtheta = squeeze( dmudthetastore(:,:,params.K,:,l,k,jcfg) );

            % calculation Poisson rates
            mu = mu+params.readnoisevariance;
            mupos = double(mu>0).*mu + double(mu<0)*keps;
            weight = 1./mupos;

            % calculation Fisher matrix
            for ii = 1:numparams
                for jj = ii:numparams
                    FisherTemp(ii,jj,l,k) = sum(weight.*dmudtheta(:,:,ii).*dmudtheta(:,:,jj),'all');
                    FisherTemp(jj,ii,l,k) = FisherTemp(ii,jj,l,k);
                end
            end

        end
    end

    Fisher = sum(FisherTemp,[3 4]);
    if (cond(Fisher)^-1>keps)
        CRLBstore(:,jcfg) = sqrt(diag(inv(Fisher+keps*eye(size(Fisher)))));
    end
    rcondstore(jcfg) = rcond(Fisher);

end


end