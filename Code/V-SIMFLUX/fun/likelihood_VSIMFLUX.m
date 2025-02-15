function [logL,gradlogL,HessianlogL] = likelihood_VSIMFLUX(params,image,mu,dmudtheta,varfit)

numparams = params.numparams;
L = size(mu,4);
K = size(mu,5);



% calculation of weight factors
keps = 1e3*eps;
mupos = double(mu>0).*mu + double(mu<0)*keps;
weight = (image-mupos)./(mupos+varfit);
dweight = (image+varfit)./(mupos+varfit).^2;

% log-likelihood, gradient vector and Hessian matrix
logL = sum((image+varfit).*log(mupos+varfit)-(mupos+varfit),1:5);

gradlogLTemp=zeros(1,numparams,L*K);
for l=1:L
    for k=1:K
        gradlogLTemp(:,:,(l-1)*K+k) = permute(sum(weight(:,:,:,l,k).*dmudtheta(:,:,:,:,l,k),1:3),[5 4 3 2 1]);
    end
end
gradlogL = sum(gradlogLTemp,[1 3]);


HessianlogLTemp = zeros(numparams,numparams,L*K);
for ii = 1:numparams
    for jj = 1:numparams
        for l=1:L
            for k=1:K
                HessianlogLTemp(ii,jj,(l-1)*K+k) = sum(-dweight(:,:,:,l,k).*dmudtheta(:,:,:,ii,l,k).*dmudtheta(:,:,:,jj,l,k),1:3);
            end
        end
    end
end
HessianlogL = sum(HessianlogLTemp,3);


end