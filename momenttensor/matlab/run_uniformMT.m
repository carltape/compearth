%
% run_uniformMT.m
%
% run script for uniformMT.m
%
% Carl Tape, 7/26/2015
%

clear, clc, close all

deg = 180/pi;

% choose examples (see descriptions below)
bex1 = true;       % random full moment tensor
bex2 = false;       % regular grid of uniform full moment tensors
bex3 = false;       % random double couple moment tensors
bex4 = false;       % regular grid of double couple moment tensors
bex5 = false;       % random moment tensors with fixed eigenvalues
bex6 = false;       % will generate error (intentionally)
banalog = false;    % insights into the sin^4(omega) distribution

% reference moment tensor for setting omega=0 in the distributions
% (alternatively you can choose a moment tensor that is in the set)
Mref = 1/sqrt(2)*[1 0 -1 0 0 0]';

if bex1
    % randomly generated uniform full moment tensors
    n = 1e5;
    [M,u,v,kappa,sigma,h] = uniformMT(n);
    %iref = randi(n); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
    plot_omega(omega);  
end

if bex2
    % regular grid of uniform full moment tensors
    %n = [2 3 5 4 2];
    n = [6 18 18 9 5];
    %n = [6 18 36 18 10];              % 10-degree increments
    [M,u,v,kappa,sigma,h] = uniformMT(n);
    ntotal = length(u);
end

if bex3
    % randomly generated uniform double couple moment tensors
    n = 1e5;
    gamma0 = 0; delta0 = 0;     % double couple
    [M,u,v,kappa,sigma,h] = uniformMT(n,gamma0,delta0);
    %iref = randi(n); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
    plot_omegadc(omega);  
end

if bex4
    % regular grid of uniform double couple moment tensors
    ntotal = 1e5;
    %ntotal = 46656;                 % number of points for n = [0 0 72 36 18]
    npart = round(ntotal^(1/3));
    %n = [0 0 5 4 2];
    %n = [0 0 36 18 9];              % 10-degree increments
    n = [0 0 72 36 18];             % 5-degree increments
    %n = [0 0 npart npart npart];    % equal number of points for each interval
    %n = [0 0 8 180 100];            % PROBLEMATIC!
    gamma0 = 0; delta0 = 0;     % double couple
    [M,u,v,kappa,sigma,h] = uniformMT(n,gamma0,delta0);
    ntotal = length(u);
    %iref = randi(ntotal); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
    plot_omegadc(omega); 
end

if bex5
    % randomly generated uniform moment tensors with fixed lune point
    n = 1e5;
    gamma0 = -25; delta0 = 60;     % double couple
    [M,u,v,kappa,sigma,h] = uniformMT(n,gamma0,delta0);
    %iref = randi(n); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
end

if bex6
    %will generate error, since the 0 0 indicates a fixed lune point,
    % but no lune point is specified as the 2nd and 3rd arguments of uniformMT
    n = [0 0 5 4 2];
    %n = [1 1 5 4 2];
    [M,u,v,kappa,sigma,h] = uniformMT(n);
end

% insights into the sin^4(omega) distribution
if banalog
    figure; nr=2; nc=1;

    % angular distance from a point on a circle (1-sphere) to all other points on a circle
    n = 1e6;
    theta = 2*pi*rand(n,1);
    [rx,ry] = pol2cart(theta,1);
    ipick = randi(n);
    rx0 = rx(ipick)*ones(n,1);
    ry0 = ry(ipick)*ones(n,1);
    theta = acos( rx0.*rx + ry0.*ry );
    subplot(nr,nc,1); plot_histo(theta*deg,[0:2:180],3); ylim([0 1e-2]);
    title('angular distance from any point on the circle to uniformly distributed points');
    %figure; hold on; plot(rx,ry,'.'); plot(rx(ipick),ry(ipick),'ro');

    % angular distance from a point on a sphere (2-sphere) to all other points
    n = 1e6;
    z = randomvec(-1,1,n);
    ylat = asin(z)*deg;         % delta = 90 - acos(z)*deg;
    xlon = randomvec(0,360,n);
    ipick = randi(n);
    xlon0 = xlon(ipick)*ones(n,1);
    ylat0 = ylat(ipick)*ones(n,1);
    % this would be cleaner if written as a dot product
    theta = distance(ylat0,xlon0,ylat,xlon,'degrees');

    subplot(nr,nc,2); hold on;
    [N,Nplot,centers] = plot_histo(theta,[0:2:180],3); ylim([0 1e-2]);
    x = centers;
    f = sin(x/deg);
    plot(x,f/2/deg,'r.-');
    title('distance from any point on the sphere to uniformly distributed points');
    xlabel('angular distance, deg');
    %figure; hold on; plot(xlon,ylat,'.'); plot(xlon(ipick),ylat(ipick),'ro');
    
    break
end

%==========================================================================

% Question: How many fewer points are needed to calculate a uniform
% distribution of moment tensors, in comparison to using assuming that the
% eigenvalue triples are uniformly distributed on the lune?
%
% Suppose he takes eigenvalue triples uniformly distributed on the lune and
% then takes the same number of random orientations at each of his lune points.
% And suppose he wants to achieve the same density of moment tensors that we
% have at the double couple. Then from our Eq 44b he wants density 4/pi
% everywhere on the lune. The integral of that over the lune is
% 4/pi x lune area = 4\pi x 4 pi^2/6.  The integral of our density is 1,
% so the other guy is using about 8 times as many points.
%
%disp(sprintf('savings in number of gridpoints is %.2f',8*pi/3));

%==========================================================================
% plot various moment tensor quantities

disp('plotting various quantities...');

n = length(kappa);
M0 = CMT2m0(1,M);
theta = acos(h)*deg;
[gamma,delta] = uv2lune(u,v);

% lune longitude, lune latitude, strike, dip, slip
plotMT_TT(gamma,delta,M0,kappa,theta,sigma);

break

% lune longitude and latitude within latitude bands and longitude bands
plotMT_lune(gamma,delta);
% Mij entries
Medges = [-sqrt(2):0.1:sqrt(2)];
plotMT_Mij(M,Medges);
% eigenvalues and plunge-azimuth angles of eigenvectors
[lam,U] = CMTdecom(M);
Useu = convertv(1,5,U);
Uout = U2pa(Useu,1);
plotMT_eigvec(Uout(:,1),Uout(:,2),Uout(:,3),Uout(:,4),Uout(:,5),Uout(:,6));
plotMT_lam(lam);

%==========================================================================
