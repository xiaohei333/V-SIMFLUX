function [Nstore,mustore,dmudthetastore,meritstore,numiters] = localization_N(allspots,theta0,params)
fprintf('\nStart fitting %i instances:\n',params.Ncfg); tic;

% parameter settings
Ncfg = params.Ncfg;
varfit = params.varfit;
tollim = params.tollim;
Nitermax = params.Nitermax;
numparams = params.numparams;

% pre-allocation
numiters = zeros(1,Ncfg);
mustore = zeros(params.Mx,params.My,params.K,Ncfg);
dmudthetastore = zeros(params.Mx,params.My,params.K,numparams,Ncfg);
Nstore = zeros(1,Ncfg,Nitermax+1);
meritstore = zeros(Ncfg,Nitermax+1);

flg_nat = params.flg_nat;
if flg_nat
    natPredictions = params.natPredictions;
else
    natPredictions = zeros(Ncfg,1);
end

% setup parallel loop
if params.flg_parallel
    p = gcp('nocreate');
    if isempty(p)
        parpool;
    end
    parforArg = Inf;
else
    parforArg = 0;
end

parfor jcfg = 1:Ncfg
    
    if flg_nat
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params,natPredictions(jcfg,:));
    else
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);
    end
    
    % pre-allocate
    Ntemp = zeros(1,Nitermax+1);
    merittemp = zeros(1,Nitermax+1);
    
    % initial values and max/min values
    theta = theta0(:,jcfg)';
    spots = allspots(:,:,:,jcfg);
    N = theta(4);

    %修改
    Nmin = N/10;
    Nmax = 2*N;
    Nretry = (Nmax+Nmin)/2;
    %修改

    [mu,dmudtheta] = poissonrate(params,theta,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
    [merit,grad,Hessian] = likelihood(params,spots,mu,dmudtheta,varfit);

    Ntemp(:,1) = N;
    merittemp(1) = merit;
    meritprev = merit;
    
    % start iteration loop
    iiter = 1;
    monitor = 2*tollim;
    alamda = 1;
    alamdafac = 1.5;

while ((iiter<=Nitermax) && (monitor>tollim))
        
        % check for det(H)=0 in order to avoid inversion of H
        matty = Hessian+alamda*diag(diag(Hessian));
        if (abs(det(matty))>2*eps)

        % update parameters    
        % update of fit parameters via Levenberg-Marquardt
        Ntry = Nupdate(N,Nmax,Nmin,Nretry,grad,Hessian,alamda);
        thetatry = theta;
        thetatry(4) = Ntry;

        % calculate update merit function
        [mu,dmudtheta] = poissonrate(params,thetatry,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
        [merittry,gradtry,Hessiantry] = likelihood(params,spots,mu,dmudtheta,varfit);
        dmerit = merittry-merit;

        % modify Levenberg-Marquardt parameter
          if (dmerit<0)
              alamda = alamdafac*alamda;
          else
              alamda = alamda/alamdafac;
              N = Ntry;
              merit = merittry;
              grad = gradtry;
              Hessian = Hessiantry;
              dmerit = merit-meritprev;
              monitor = abs(dmerit/merit);
              meritprev = merit;
              Nretry = N;
          end

        else
            alamda = alamdafac*alamda;
        end
     
        % store values and update counter
        Ntemp(:,iiter+1) = N;
        merittemp(iiter+1) = merit;
        
        iiter = iiter+1; % update counter

end

    % store values
    numiters(jcfg) = iiter-1;
    
    for jiter = iiter+1:Nitermax+1
        merittemp(jiter) = merit;
        Ntemp(:,jiter) = N;
    end
    
    mustore(:,:,:,jcfg) = mu;
    dmudthetastore(:,:,:,:,jcfg) = dmudtheta;
    Nstore(:,jcfg,:) = Ntemp;
    meritstore(jcfg,:) = merittemp;
    
    % print update
    if rem(jcfg,round(Ncfg/10)) == 0
        fprintf('fitting instance # %i...\n',jcfg)
    end

end



% add offset to merit function
meritoffset = meritoffsetcalc(allspots,params.varfit);
for jiter = 1:Nitermax+1
    meritstore(:,jiter) = meritstore(:,jiter)+meritoffset;
end

% print run time
fprintf(['\nMLE fit routine (spot/second): ' num2str(toc,3) 's (' num2str(params.Ncfg/toc,5) ')\n'])


end



function Nnew = Nupdate(Nold,Nmax,Nmin,Nretry,grad,Hessian,alamda)

        Bmat = Hessian(4,4)+alamda*diag(diag(Hessian(4,4)));
        dN = -Bmat\grad(4)';
        Nnew = Nold+dN';
        % enforce physical boundaries in parameter space.       
        if ((Nnew>Nmax)||(Nnew<Nmin))
                    Nnew = Nretry;
        end

end
