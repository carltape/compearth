%
% forward_epicenter.m
%
% This forward problem describes the straight-line ray path travel time
% from a source with an unknown origin time to a receiver in a
% homogeneous medium with unknown velocity V.
%
% NOTE: variables such as minitial, mtarget, and eobs may be over-written
% within optimization.m when the user specifies multiple runs.
%
% NOTE: If there are constants in the Matlab in-line functions, then these
% constants are assigned at the time the function is INITIALIZED, not at
% the time the function is called.
%
% CALLED BY : optimization.m
% CALLS     : plot_epicenters.m
%
% Carl Tape, 3/18/2010
%

%==========================================================================
% FORWARD PROBLEM    

ndata = 12;
nparm = 4;

% labels for parameters
%mlabs = {'xs','ys','ts','v = ln(V/V0)'};
mlabs = {'xs','ys','ts','v'};
ulabs = {'km','km','s','none'};

% velocity for prior model
% v is the unitless logarithmic velocity, v = ln(V/V0)
% V = V0*exp(v) has units km/s
V = 5;      % km/s
V0 = 1;     % arbitrary scale factor (V is defined below)

% travel time computation (homogeneous velocity; straight ray paths)
d2 = @(x,y,xr,yr)     ( (xr-x)^2 + (yr-y)^2 );
d1 = @(x,y,xr,yr)     ( sqrt(d2(x,y,xr,yr)) );
tt = @(x,y,t,v,xr,yr) ( t + d1(x,y,xr,yr)/(V0*exp(v)) );

% receiver locations for uniform grid
xrecmin = 10; xrecmax = 80;
yrecmin = 20; yrecmax = 90;
xvec = linspace(xrecmin,xrecmax,4);
yvec = linspace(yrecmin,yrecmax,3);
[X,Y] = meshgrid(xvec,yvec);
xrec = X(:);
yrec = Y(:);

% computation of predictions
d = @(m) ([   tt(m(1),m(2),m(3),m(4),xrec(1),yrec(1))
              tt(m(1),m(2),m(3),m(4),xrec(2),yrec(2))
              tt(m(1),m(2),m(3),m(4),xrec(3),yrec(3))
              tt(m(1),m(2),m(3),m(4),xrec(4),yrec(4))
              tt(m(1),m(2),m(3),m(4),xrec(5),yrec(5))
              tt(m(1),m(2),m(3),m(4),xrec(6),yrec(6))
              tt(m(1),m(2),m(3),m(4),xrec(7),yrec(7))
              tt(m(1),m(2),m(3),m(4),xrec(8),yrec(8))
              tt(m(1),m(2),m(3),m(4),xrec(9),yrec(9))
              tt(m(1),m(2),m(3),m(4),xrec(10),yrec(10))
              tt(m(1),m(2),m(3),m(4),xrec(11),yrec(11))
              tt(m(1),m(2),m(3),m(4),xrec(12),yrec(12))
           ]);

% N x M matrix of partial derivatives (differentiate tt with respect to each parameter)
% evaluated at model m = (m(1),m(2),m(3),m(4))
G = @(m) ([
    -(d2(m(1),m(2),xrec(1),yrec(1)))^(-1/2)*(xrec(1)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(1),yrec(1)))^(-1/2)*(yrec(1)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(1),yrec(1))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(2),yrec(2)))^(-1/2)*(xrec(2)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(2),yrec(2)))^(-1/2)*(yrec(2)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(2),yrec(2))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(3),yrec(3)))^(-1/2)*(xrec(3)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(3),yrec(3)))^(-1/2)*(yrec(3)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(3),yrec(3))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(4),yrec(4)))^(-1/2)*(xrec(4)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(4),yrec(4)))^(-1/2)*(yrec(4)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(4),yrec(4))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(5),yrec(5)))^(-1/2)*(xrec(5)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(5),yrec(5)))^(-1/2)*(yrec(5)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(5),yrec(5))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(6),yrec(6)))^(-1/2)*(xrec(6)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(6),yrec(6)))^(-1/2)*(yrec(6)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(6),yrec(6))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(7),yrec(7)))^(-1/2)*(xrec(7)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(7),yrec(7)))^(-1/2)*(yrec(7)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(7),yrec(7))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(8),yrec(8)))^(-1/2)*(xrec(8)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(8),yrec(8)))^(-1/2)*(yrec(8)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(8),yrec(8))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(9),yrec(9)))^(-1/2)*(xrec(9)-m(1))/(V0*exp(m(4)))       -(d2(m(1),m(2),xrec(9),yrec(9)))^(-1/2)*(yrec(9)-m(2))/(V0*exp(m(4)))       1  -d1(m(1),m(2),xrec(9),yrec(9))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(10),yrec(10)))^(-1/2)*(xrec(10)-m(1))/(V0*exp(m(4)))    -(d2(m(1),m(2),xrec(10),yrec(10)))^(-1/2)*(yrec(10)-m(2))/(V0*exp(m(4)))    1  -d1(m(1),m(2),xrec(10),yrec(10))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(11),yrec(11)))^(-1/2)*(xrec(11)-m(1))/(V0*exp(m(4)))    -(d2(m(1),m(2),xrec(11),yrec(11)))^(-1/2)*(yrec(11)-m(2))/(V0*exp(m(4)))    1  -d1(m(1),m(2),xrec(11),yrec(11))/(V0*exp(m(4)))
    -(d2(m(1),m(2),xrec(12),yrec(12)))^(-1/2)*(xrec(12)-m(1))/(V0*exp(m(4)))    -(d2(m(1),m(2),xrec(12),yrec(12)))^(-1/2)*(yrec(12)-m(2))/(V0*exp(m(4)))    1  -d1(m(1),m(2),xrec(12),yrec(12))/(V0*exp(m(4)))
    ]);

% M x M matrix of second partial derivatives (only used in full Newton method)
% note: this contains the measurement index ii
G2 = @(m,ii)  ([
   (d1(m(1),m(2),xrec(ii),yrec(ii)))^-3*(yrec(ii)-m(2))^2/(V0*exp(m(4)))               -(d1(m(1),m(2),xrec(ii),yrec(ii)))^-3*(xrec(ii)-m(1))*(yrec(ii)-m(2))/(V0*exp(m(4)))     0  (d1(m(1),m(2),xrec(ii),yrec(ii)))^-1*(xrec(ii)-m(1))/(V0*exp(m(4)))
  -(d1(m(1),m(2),xrec(ii),yrec(ii)))^-3*(xrec(ii)-m(1))*(yrec(ii)-m(2))/(V0*exp(m(4)))   (d1(m(1),m(2),xrec(ii),yrec(ii)))^-3*(xrec(ii)-m(1))^2/(V0*exp(m(4)))                  0  (d1(m(1),m(2),xrec(ii),yrec(ii)))^-1*(yrec(ii)-m(2))/(V0*exp(m(4)))
    0                                                                                0                                                                                   0   0
   (d1(m(1),m(2),xrec(ii),yrec(ii)))^-1*(xrec(ii)-m(1))/(V0*exp(m(4)))                  (d1(m(1),m(2),xrec(ii),yrec(ii)))^-1*(yrec(ii)-m(2))/(V0*exp(m(4)))                    0   d1(m(1),m(2),xrec(ii),yrec(ii))/(V0*exp(m(4)))
]);

%--------------------------------------------------------------------------
% RANDOM VECTORS (for sampling covariance matrices)

% Gaussian random vectors, each with mean = 0 and standard deviation = 1
randn_vecs_m = randn(nparm,nsamples);   % model
randn_vecs_d = randn(ndata,nsamples);   % data

%--------------------------------------------------------------------------
% PRIOR MODEL (MEAN MODEL) : ts, xs, ys, v

% prior model
% note: V and V0 are defined above
mprior = [ 35 45 16 log(V/V0) ]';

% prior model covariance matrix (assumed to be diagonal)
sigma_prior = [10 10 0.5 0.2]';         % standard deviations
cprior0     = diag( sigma_prior.^2 );   % diagonal covariance matrix
if inormalization==1
    Cmfac = nparm;
else
    Cmfac = 1;
end
cprior   = Cmfac * cprior0;             % WITH NORMALIZATION FACTOR
icprior  = inv(cprior);                 % WITH NORMALIZATION FACTOR
icprior0 = inv(cprior0);
Lprior   = chol(cprior0,'lower')';      % square-root (lower triangular)

% sample the prior model distribution using the square-root UNNORMALIZED covariance matrix
cov_samples_m  = Lprior * randn_vecs_m;
mprior_samples = repmat(mprior,1,nsamples) + cov_samples_m;

% compute the norm of each model sample using the inverse NORMALIZED covariance matrix
norm2_mprior = zeros(nsamples,1);
for xx=1:nsamples
    dm = mprior_samples(:,xx) - mprior;
    norm2_mprior(xx) = dm' * icprior * dm;
end
%figure; plot(norm2_mprior,'.')

%--------------------------------------------------------------------------
% INITIAL MODEL

% minitial is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
if irandom_initial_model == 1
    minitial = mprior_samples(:,1);     % first sample (random)
else
    minitial = [                        % fixed
        46.5236
        40.1182
        15.3890
        1.7748
    ];
end

%--------------------------------------------------------------------------
% TARGET MODEL

% mtarget is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
if irandom_target_model == 1  
    mtarget = mprior_samples(:,end);    % last sample (random)
else
    mtarget = [                         % fixed  
           21.2922
           46.2974
           16.1314
            2.0903
                ];
end

%--------------------------------------------------------------------------

% TARGET DATA
dtarget = d(mtarget);

% data covariance matrix (assumed to be diagonal)
tsigma = 0.5;                       % uncertainty in arrival time measurement, seconds
sigma_obs = tsigma * ones(ndata,1); % standard deviations
cobs0     = diag( sigma_obs.^2 );   % diagonal covariance matrix
if inormalization==1
    Cdfac = ndata;
else
    Cdfac = 1;
end
cobs   = Cdfac * cobs0;             % WITH NORMALIZATION FACTOR
icobs  = inv(cobs);                 % WITH NORMALIZATION FACTOR
icobs0 = inv(cobs0);
Lcobs  = chol(cobs0,'lower')';      % square-root (lower triangular)

% sample the data distribution using the square-root UNNORMALIZED covariance matrix
cov_samples_d = Lcobs * randn_vecs_d;
dobs_samples  = repmat(dtarget,1,nsamples) + cov_samples_d;

% compute the norm of each data sample using the inverse NORMALIZED covariance matrix
norm2_dobs = zeros(nsamples,1);
for xx=1:nsamples
    dd = dobs_samples(:,xx) - dtarget;
    norm2_dobs(xx) = dd' * icobs * dd;
end
%figure; plot(norm2_dobs,'.')

% Pick the uncertainties for the target data by simply choosing a sample
%   from the realizations of the data covariance.
% NOTE: eobs is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
switch idata_errors
    case 0
        eobs = zeros(ndata,1);          % no errors
    case 1
        eobs = cov_samples_d(:,1);      % first sample (random)
    case 2
        eobs = [                        % fixed
            -0.8689
            -0.4666
            -0.0516
            0.2411
            0.2825
            -0.2301
            0.1977
            0.3291
            1.0063
            0.5674
            0.1348
            0.4603
        ];
end

% "true" observations (includes added errors)
dobs = dtarget + eobs;

%--------------------------------------------------------------------------
% PLOTS

axepi = [0 100 0 100];

% source-receiver geometry with ray paths
plot_epicenters([],mprior,minitial,mtarget,{xrec,yrec,1,axepi});
if iprint==1, print(gcf,'-depsc',sprintf('%srcrec_rays_f%i',pdir,iforward)); end

% source-receiver geometry with ray paths
plot_epicenters(mprior_samples,mprior,minitial,mtarget,{xrec,yrec,1,axepi});
if iprint==1, print(gcf,'-depsc',sprintf('%smprior_rays_f%i',pdir,iforward)); end

% with prior samples (no ray paths)
opts = {xrec,yrec,0,axepi};
plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts);
if iprint==1, print(gcf,'-depsc',sprintf('%smprior_f%i',pdir,iforward)); end

%==========================================================================
