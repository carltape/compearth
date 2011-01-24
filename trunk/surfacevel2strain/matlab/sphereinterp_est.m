%
% function sphereinterp_est
% Carl Tape and Pablo Muse, 04-Jan-2011
%
% Estimate a smooth field on the sphere from discrete points using
% spherical wavelets. This is a stripped-down version of the procedure
% presented in surfacevel2strain (Tape et al. 2009, GJI).
%
% INPUT
%   dlon    longitude of data points
%   dlat    latitude of data points
%   d    scalar value of data points
%   dsig    estimated uncertainties associated with each scalar value
%   ax0     lon-lat box contained desired region
%   pparm   plotting parameters
%        1. nx, number of gridpoints in x direction
%        2. ulabel, label for the observations (e.g., 'zdep, km')
%        3. polyxlon, polygon containing plotting points (optional)
%        4. polyylat, polygon containing plotting points (optional)
%
% calls dogsph_vals.m, ridge_carl.m
% called by sphereinterp.m
%

function [dest,dest_plot,lam0,dlon_plot,dlat_plot,na,nb] ...
    = sphereinterp_est(spline_tot,dlon,dlat,d,dsig,ax0,rparm,pparm)

disp('------------------------------------------------------------');
disp('entering sphereinterp_est.m to estimate smooth scalar field on the sphere');

earthr = 6371*1e3;      % earth radius (m)
deg = 180/pi;
msize = 6^2;
ngrid = length(spline_tot);

dlab = 'zdep (km)';  % USER PARAMETER

nlam = rparm{1};
ilampick = rparm{2};  % =1 (iL), =2 (iOCV), =3 (iGCV)

numx_plot = pparm{1};
ulabel = pparm{2};
ipoly = 0;
if length(pparm) > 2
    polylon = pparm{3};
    polylat = pparm{4};
    disp(sprintf('input polygon has %i points',length(polylon)));
    figure; plot(polylon,polylat,'.-');
    ipoly = 1;
else
    disp('no input polygon provided --> full plotting grid will be used');
end

lonmin = ax0(1); lonmax = ax0(2);
latmin = ax0(3); latmax = ax0(4);
ndata = length(dlon);

disp('choice of regularization parameter:');
if any([1 2 3]==ilampick)
    if ilampick==1, streg = 'L-curve'; end
    if ilampick==2, streg = 'ordinary cross validation'; end
    if ilampick==3, streg = 'generalized cross validation'; end
elseif ilampick < 0
    streg = 'manually picked parameter';
    if or(abs(ilampick) < 1, abs(ilampick) > nlam)
        error(sprintf('ilampick (%i) must be 1 <= ilam <= %i',abs(ilampick),nlam));
    end
else
    error(sprintf('ilampick (%i) must be 1 (iL), 2 (iOCV), or 3 (iGCV)',ilampick));
end
disp(streg);

%-------------------------------
% USER PARAMETERS

slabel = 'sphereinterp_est.m';

% min number of gridpoints for a particular order
%nmin = 1;

% threshold wavelets based on data
%ntrsh = 3;       % KEY: number of evaluations >= qtrsh

% =1 to use weights, =0 to ignore weights
if isempty(dsig)
    disp('no input uncertainties provided');
    icov = 0;
else
    disp('input uncertainties will be used within inversion');
    icov = 1;
end

% q for the "secular field" -- which combines the estimates from qmin to
% qsec for the multiscale analysis (NOT VERY RELEVANT HERE)
%qsec = round(mean([qmin qmax]));

% vector of regularization (damping) parameters
% NOTE: in many cases, you MUST have a non-zero damping parameter
%nlam = 40;
if icov==0          % unweighted
    minlampwr = -8; maxlampwr = 2;
else                % weighted
    minlampwr = -3; maxlampwr = 6;
end
lampwr = linspace(minlampwr,maxlampwr,nlam);
lamvec = 10.^lampwr;

%========================================================
% PROPERTIES OF THE GRIDS AND THE SPHERICAL WAVELETS


%========================================================
% DAMPED LSQ TO ESTIMATE THE VELOCITY FIELD

% REGULARIZATION MATRIX (MODEL COVARIANCE MATRIX)
% diagonal of regularization matrix (or model covariance matrix)
scales = spline_tot(:,3);   % CHT 8/3/08, changed from spline_tot(:,3)-1
pow2scales = 2.^scales;                 % is it 2^p or 4^p ?
Dmat = diag(pow2scales.^2);
Dhalfinv = diag(1 ./ pow2scales );      % regularization by gradient-model norm
%Dhalfinv = eye(ngrid);                  % regularization by model-norm

%--------------------------------------------------------
% CONSTRUCT DESIGN MATRIX for components of velocity field

disp('  '); disp('Constructing the design matrix...');

% get the 'base' design matrix
[G, Gdph, Gdth] = dogsph_vals_mat(spline_tot, dlon, dlat, {3});

% number of regularization parameters
lam0 = NaN;       % lam0(1) = NaN if ndim = 2

%-------------------------------

stlams = [' lam = ' num2str(sprintf('%.2f',lamvec(1)))  ...
    ' to ' num2str(sprintf('%.2f',lamvec(end))) '  (' num2str(nlam) ' solutions)'];

% data vector and weighting vector (looping order is r-theta-phi)
% DATA COVARIANCE MATRIX
d = d; wu = 1 ./ dsig.^2; Wvec = wu;

if icov == 0
    wu = ones(ndata,1);
    ws = ones(ndata,1);
    we = ones(ndata,1);
    Wvec = ones(ndata,1);
end

disp(' creating the L-curve...');
trms = zeros(nlam,1);
mss = zeros(nlam,1);

Gvec = zeros(nlam,1);
rss0 = 0; mss0 = 0; G0 = 0;
Whalf = diag( sqrt(Wvec) );     % Weisberg, p. 97
[f_h_prime, rss, mss, Gvec, Fvec, dof, kap, iL, iGCV, iOCV] = ...
    ridge_carl(Whalf*d, Whalf*G*Dhalfinv, lamvec);

% (un-)transform model vector
f_h = zeros(ngrid,nlam);
for ik = 1:nlam
    f_h(:,ik) = Dhalfinv * f_h_prime(:,ik);
end

% KEY: select on the basis of the OCV curve, GCV curve, or L curve
if ilampick==1
    ilam = iL;
elseif ilampick==2
    ilam = iOCV;
elseif ilampick==3
    ilam = iGCV;
else
    ilam = abs(ilampick);
end
lam0 = lamvec(ilam);

disp('  ');
disp('Pick the regularization parameter:');
disp(sprintf('L-curve lambda = %.3e (index %i)',lamvec(iL),iL));
disp(sprintf('    OCV lambda = %.3e (index %i)',lamvec(iOCV),iOCV));
disp(sprintf('    GCV lambda = %.3e (index %i)',lamvec(iGCV),iGCV));
disp(sprintf('your pick lam0 = %.3e (index %i)',lam0,ilam));

%========================================================

disp(' computing the model vector...');

%--------------------------------
% COMPUTE MODEL VECTOR

fu = zeros(ngrid,1);

% least squares solution
Cm_u = inv(G'*diag(wu)*G + lam0^2*Dmat);   % m x m
fu = Cm_u*G'*diag(wu)*d;                   % m x 1
dest = G*fu;                               % n x 1
d_res = d - dest;                          % n x 1
Cd_u_diag = diag(G*Cm_u*G');               % n x 1
dsig_post = sqrt(Cd_u_diag);               % n x 1

figure; nr=2; nc=2;
resmax = max(abs(d_res));

subplot(2,1,1); hold on;
plot( d, dest, '.');
%plot(max(abs(d))*[-1 1],max(abs(d))*[-1 1],'r--');
plot([0 60],[0 60],'r--');
xlabel(['OBSERVED ' ulabel]);
ylabel(['ESTIMATED ' ulabel]);
%title(['cor(d-obs, d-est) = ' num2str(corr(d,dest))]);  %
%STATISTICS TOOLBOX -- corr
axis equal, grid on;

subplot(nr,nc,3); hold on;
plot(d_res, '.');
xlabel('Observation number ');
ylabel(['RESIDUALS (d - dest), ' ulabel]);
title(['median(abs(res)) = ' num2str(median(abs(d_res))) ' km']);
grid on;

edges = linspace(-resmax,resmax,15);
subplot(nr,nc,4); plot_histo(d_res,edges);
xlabel(['RESIDUALS (d - dest), ' ulabel]);

orient tall, wysiwyg, fontsize(9)

% colored scatterplot map of residuals
figure; hold on;
scatter(dlon,dlat,msize,d_res,'filled');
caxis([-1 1]*0.5*max(abs(d_res))); axis(ax0); grid on;
colorbar; xlabel('Longitude'); ylabel('Latitude');
title(['RESIDUALS (d - dest), ' ulabel]);

%========================================================
% UNIFORM MESH FOR PLOTTING SCALAR FIELDS (from irregular data)
% (compute new design matrices with more rows but same number of columns)

% design matrix for UNIFORM plotting grid of scalar quantities
%numx_plot = input([' Plotting grid, number of points in x direction (try 200): ']);

[dlon_plot,dlat_plot,numy_plot,na,nb] = gridvec(lonmin,lonmax,numx_plot,latmin,latmax);

% only plot points inside the designated region
% NOTE: this could be done in UTM coordinates for improved accuracy
if ipoly==1
    in_plot = inpolygon(dlon_plot,dlat_plot,polylon,polylat);
    dlon_plot = dlon_plot(in_plot);
    dlat_plot = dlat_plot(in_plot);
end
nplot = length(dlon_plot);

figure; plot(dlon_plot,dlat_plot,'.'); axis(ax0);
title(sprintf('%i x %i = %i plotting points in uniform grid',numx_plot,numy_plot,nplot));

disp('  '); disp(' computing values at the plotting points...');

% fill each column of G with a basis function evaluated at all the datapoints
% NOTE: This is super slow, but it avoids having to compute Gplot (nplot x
%       ngrid) in full, as done in surfacevel2strain.m
%       A more efficient algorithm is needed.
dest_plot = zeros(nplot,1);
tic
for ii=1:nplot          % loop over data points
    if mod(ii,100)==0, disp(sprintf('%i out of %i',ii,nplot)); end
    Grow = zeros(1,ngrid);
    for jj=1:ngrid      % loop over basis functions
        Grow(jj) = dogsph_vals(spline_tot(jj,1), spline_tot(jj,2), spline_tot(jj,3), dlon_plot(ii), dlat_plot(ii), {1});
    end
    dest_plot(ii) = Grow * fu;
end
toc

% estimated field (on plotting grid)
figure; hold on;
if ipoly==1
    scatter(dlon_plot,dlat_plot,msize,dest_plot,'filled');
    %scatter(dlon_plot,dlat_plot,msize,'ko');   % black outline
else
    X = reshape(dlon_plot,na,nb);
    Y = reshape(dlat_plot,na,nb);
    Z = reshape(dest_plot,na,nb);
    pcolor(X,Y,Z); shading interp;
    scatter(dlon,dlat,msize,d,'filled');
    scatter(dlon,dlat,msize,'ko');
end
colorbar; axis(ax0);
xlabel('Longitude'); ylabel('Latitude');
title(['Estimated field (plotting grid), ' ulabel]);

%========================================================
