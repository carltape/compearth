%
% forward_epicenter_crescent.m
%
% Epicenter problem in Tarantola (2005), Problem 7-1.
%
% NOTE: If there are constants in the Matlab in-line functions, then these
% constants are assigned at the time the function is INITIALIZED, not at
% the time the function is called.
%
% Carl Tape, 23-Feb-2013
%
% CALLED BY : optimization.m
% CALLS     : plot_epicenters.m
%

%=========================================
% FORWARD PROBLEM    

ndata = 6;
nparm = 2;

% labels for parameters
mlabs = {'xs','ys'};
ulabs = {'km','km'};

% station locations, in km
xrec = [3 3 4 4 5 5];
yrec = [15 16 15 16 15 16];

% travel time computation (homogeneous velocity; straight ray paths)
V = 5;      % km/s
d2 = @(x,y,xr,yr) ( (xr-x)^2 + (yr-y)^2 );
d1 = @(x,y,xr,yr) ( sqrt(d2(x,y,xr,yr)) );
tt = @(x,y,xr,yr) ( d1(x,y,xr,yr)/V );

% computation of predictions
d = @(m) ([   tt(m(1),m(2),xrec(1),yrec(1))
              tt(m(1),m(2),xrec(2),yrec(2))
              tt(m(1),m(2),xrec(3),yrec(3))
              tt(m(1),m(2),xrec(4),yrec(4))
              tt(m(1),m(2),xrec(5),yrec(5))
              tt(m(1),m(2),xrec(6),yrec(6))
          ]);

% N x M matrix of partial derivatives (differentiate tt with respect to each parameter)
% evaluated at model m = (t,x,y,v)
G = @(m) ([
    -(d2(m(1),m(2),xrec(1),yrec(1)))^(-1/2)*(xrec(1)-m(1))/V   -(d2(m(1),m(2),xrec(1),yrec(1)))^(-1/2)*(yrec(1)-m(2))/V
    -(d2(m(1),m(2),xrec(2),yrec(2)))^(-1/2)*(xrec(2)-m(1))/V   -(d2(m(1),m(2),xrec(2),yrec(2)))^(-1/2)*(yrec(2)-m(2))/V
    -(d2(m(1),m(2),xrec(3),yrec(3)))^(-1/2)*(xrec(3)-m(1))/V   -(d2(m(1),m(2),xrec(3),yrec(3)))^(-1/2)*(yrec(3)-m(2))/V
    -(d2(m(1),m(2),xrec(4),yrec(4)))^(-1/2)*(xrec(4)-m(1))/V   -(d2(m(1),m(2),xrec(4),yrec(4)))^(-1/2)*(yrec(4)-m(2))/V
    -(d2(m(1),m(2),xrec(5),yrec(5)))^(-1/2)*(xrec(5)-m(1))/V   -(d2(m(1),m(2),xrec(5),yrec(5)))^(-1/2)*(yrec(5)-m(2))/V
    -(d2(m(1),m(2),xrec(6),yrec(6)))^(-1/2)*(xrec(6)-m(1))/V   -(d2(m(1),m(2),xrec(6),yrec(6)))^(-1/2)*(yrec(6)-m(2))/V
    ]);

% M x M matrix of second partial derivatives (only used in full Newton method)
% note: this contains the measurement index i
G2 = @(m,i)  ([
   (d1(m(1),m(2),xrec(i),yrec(i)))^-3*(yrec(i)-m(2))^2/V                -(d1(m(1),m(2),xrec(i),yrec(i)))^-3*(xrec(i)-m(1))*(yrec(i)-m(2))/V    
  -(d1(m(1),m(2),xrec(i),yrec(i)))^-3*(xrec(i)-m(1))*(yrec(i)-m(2))/V   (d1(m(1),m(2),xrec(i),yrec(i)))^-3*(xrec(i)-m(1))^2/V                  
]);

%---------------------------------------------
% RANDOM VECTORS (for sampling covariance matrices)

% Gaussian random vectors, each with mean = 0 and standard deviation = 1
randn_vecs_m = randn(nparm,nsamples);   % model
randn_vecs_d = randn(ndata,nsamples);   % data
%for ii=1:nsamples, randn_vecs_m(:,ii) = randn(nparm,1); end  % model
%for ii=1:nsamples, randn_vecs_d(:,ii) = randn(ndata,1); end  % data

%---------------------------------------------
% PRIOR MODEL (MEAN MODEL) : ts, xs, ys, v
               
% range of model space (Tarantola Figure 7.1)
xmin = 0; xmax = 22;
ymin = -2; ymax = 30;
xcen = (xmax+xmin)/2;
ycen = (ymin+ymax)/2;

% prior model
mprior = [ xcen ycen ]';

% prior model covariance matrix (diagonal)
sigma_prior = [10 10]';                 % standard deviations
cprior0     = diag( sigma_prior.^2 );   % diagonal covariance matrix
if inormalization==1
    Cmfac = nparm;
else
    Cmfac = 1;
end
cprior      = Cmfac * cprior0;          % WITH NORMALIZATION FACTOR
icprior     = inv(cprior);              % WITH NORMALIZATION FACTOR
icprior0    = inv(cprior0);
Lprior      = chol(cprior0,'lower')';   % square-root (lower triangular)

% sample the prior model distribution using the square-root UNNORMALIZED covariance matrix
for ii=1:nsamples, randn_vecs_m(:,ii) = randn(nparm,1); end
cov_samples_m  = Lprior * randn_vecs_m;
mprior_samples = repmat(mprior,1,nsamples) + cov_samples_m;

% compute the norm of each model sample using the inverse NORMALIZED covariance matrix
norm2_mprior = zeros(nsamples,1);
for ii=1:nsamples
    dm = mprior_samples(:,ii) - mprior;
    norm2_mprior(ii) = dm' * icprior * dm;
end
%figure; plot(norm2_mprior,'.')

%---------------------------------------------
% INITIAL MODEL

% minitial is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
if irandom_initial_model == 1
    minitial = mprior_samples(:,1);      % first sample
else
    minitial = [
        15
        20
    ];
end

%---------------------------------------------
% TARGET MODEL

% mtarget is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
if irandom_target_model == 1  
    mtarget = mprior_samples(:,end);      % last sample
else
    % it looks like he picked (15,5) as the true model,
    % but he did NOT add errors consistent with tsigma=0.1 s 
    %tobs = tm(xtar,ytar) + tsigma*randn(N,1);
    mtarget = [
           15
            5
                ];
end

%---------------------------------------------

% TARGET DATA
dtarget = d(mtarget);

% data covariance matrix
tsigma = 0.1;                       % uncertainty in travel time measurement, seconds
sigma_obs = tsigma * ones(ndata,1); % standard deviations
cobs0     = diag( sigma_obs.^2 );   % diagonal covariance matrix
if inormalization==1
    Cdfac = ndata;
else
    Cdfac = 1;
end
cobs      = Cdfac * cobs0;          % WITH NORMALIZATION FACTOR
icobs     = inv(cobs);              % WITH NORMALIZATION FACTOR
icobs0    = inv(cobs0);
Lcobs     = chol(cobs0,'lower')';   % square-root (lower triangular)

% sample the data distribution using the square-root UNNORMALIZED covariance matrix
for ii=1:nsamples, randn_vecs_d(:,ii) = randn(ndata,1); end
cov_samples_d = Lcobs * randn_vecs_d;
dobs_samples  = repmat(dtarget,1,nsamples) + cov_samples_d;

% compute the norm of each data sample using the inverse NORMALIZED covariance matrix
norm2_dobs = zeros(nsamples,1);
for ii=1:nsamples
    dd = dobs_samples(:,ii) - dtarget;
    norm2_dobs(ii) = dd' * icobs * dd;
end
%figure; plot(norm2_dobs,'.')

% Pick the uncertainties for the target data by simply choosing a sample
%   from the realizations of the data covariance.
% NOTE: eobs is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
switch idata_errors
    case 0
        eobs = zeros(ndata,1);          % no errors
    case 1
        eobs = cov_samples_d(:,1);      % first sample
    case 2
        %tobs = [3.12 3.26 2.98 3.12 2.84 1]';       % outlier for 6th station
        tobs = [3.12 3.26 2.98 3.12 2.84 2.98]';    % fixed
        eobs = tobs - dtarget;
        % note: these are MUCH less than tsigma = 0.1 and suggest that the
        %       added errors are inconsistent with the data covariance matrix
        %eobs = [
        %   -0.0041
        %    0.0042
        %    0.0068
        %    0.0087
        %    0.0116
        %    0.0068 ];
end

% "true" observations (includes added errors)
dobs = dtarget + eobs;

%---------------------------------------------
% PLOTS

axepi = [xmin xmax ymin ymax];

% source-receiver geometry with ray paths
plot_epicenters([],mprior,minitial,mtarget,{xrec,yrec,1,axepi});
if iprint==1, print(gcf,'-depsc',sprintf('%ssrcrec_rays_f%i',pdir,iforward)); end

% source-receiver geometry with ray paths
plot_epicenters(mprior_samples,mprior,minitial,mtarget,{xrec,yrec,1,axepi});
if iprint==1, print(gcf,'-depsc',sprintf('%smprior_rays_f%i',pdir,iforward)); end

% with prior samples (no ray paths)
opts = {xrec,yrec,0,axepi};
plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts);
if iprint==1, print(gcf,'-depsc',sprintf('%smprior_f%i',pdir,iforward)); end

%========================================================================
