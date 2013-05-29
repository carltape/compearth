function [gamma, epsilon, eCLVD, eDC, trM, eISO, eDEV] = CMT2epsilon(M,iepsilon,isort)
% CMT2EPSILON converts a moment tensor to CLVD-related parameters
%
% This function converts from a (CMT) moment tensor to epsilon, which
% represents the magnitude of the non-double-couple component in a moment
% tensor.
%
% INPUT
%   M           6 x n set of input moment tensors
%   iepsilon    OPTIONAL: choice of equation for epsilon
%   isort       OPTIONAL: choice of sorting of eigenvales (CMTdecom.m)
%                   1: lam1 >= lam2 >= lam3
%                   2: lam1 <= lam2 <= lam3
%                   3: | lam1 | >= | lam2 | >= | lam3 |
%                   4: | lam1 | <= | lam2 | <= | lam3 |
%
% Equation for epsilon (see Latex notes):
%    Dziewonski et al. (1981), p. 2837
%    Jost and Herrmann (1989), Eq. 38
%
% moment tensor M = [Mrr Mtt Mpp Mrt Mrp Mtp]
%
% calls Mdim.m, CMTdecom_iso.m, CMTdecom.m
%
% Carl Tape, 01-April-2011
%

% default choices
if nargin==1
    iepsilon = 1;   % choice of epsilon expression (see below)
    isort = 1;      % choice of eigenvalue sorting (CMTdecom.m)
end

% make sure M is 6 x n
[M,n] = Mdim(M);

Mrr = M(1,:); Mtt = M(2,:); Mpp = M(3,:);
Mrt = M(4,:); Mrp = M(5,:); Mtp = M(6,:);

% decomposition into isotropic and deviatoric parts
[Miso,Mdev,trM] = CMTdecom_iso(M);
lamiso = Miso(1:3,:);

% decomposition of DEVIATORIC moment tensor into eigenbasis
% NOTE: ordering of eigenvalues is important
[lamdev,~] = CMTdecom(Mdev,isort);

% epsilon: a measure of CLVD relative magnitude
% NOTE: each case takes into account a different sorting of eigenvalues
if iepsilon == 1
    % traditional choice ("double projection") -- range from [-0.5 to 0.5]
    switch isort
        case 1, epsilon = -lamdev(2,:) ./ max(abs(lamdev([1 3],:)));
        case 2, epsilon = -lamdev(2,:) ./ max(abs(lamdev([1 3],:)));
        case 3, epsilon = -lamdev(3,:) ./ max(abs(lamdev([1 2],:))); 
        case 4, epsilon = -lamdev(1,:) ./ max(abs(lamdev([2 3],:))); 
    end
    
elseif iepsilon == 2
    % traditional choice ("double projection") -- range from [0 to 0.5]
    switch isort
        case 1, epsilon = abs(lamdev(2,:)) ./ max(abs(lamdev([1 3],:)));
        case 2, epsilon = abs(lamdev(2,:)) ./ max(abs(lamdev([1 3],:)));
        case 3, epsilon = abs( lamdev(3,:) ./ lamdev(1,:) ); 
        case 4, epsilon = abs( lamdev(1,:) ./ lamdev(3,:) ); 
    end
    
elseif iepsilon == 3  
    % "single projection choice" -- range from [-1.0 to 0.5]
	switch isort
        case 1, epsilon = -lamdev(2,:) ./ lamdev(3,:);
        case 2, epsilon = -lamdev(2,:) ./ lamdev(1,:);
        case 3, epsilon = -lamdev(3,:) ./ min(lamdev(1:2,:));
        case 4, epsilon = -lamdev(1,:) ./ min(lamdev(2:3,:)); 
    end
    
elseif iepsilon == 4    
    % Walt's choice -- such that gamma, epsilon, lam2 all have the same sign
    % note: same as iepsilon=1 but with flipped sign
    switch isort
        case 1, epsilon = lamdev(2,:) ./ max(abs(lamdev([1 3],:)));
        case 2, epsilon = lamdev(2,:) ./ max(abs(lamdev([1 3],:)));
        case 3, epsilon = lamdev(3,:) ./ max(abs(lamdev([1 2],:))); 
        case 4, epsilon = lamdev(1,:) ./ max(abs(lamdev([2 3],:))); 
    end   
    
else
    error(sprintf('iepsilon (%i) must be 1,2, or 3',iepsilon));
end
epsilon = epsilon(:);

% compute gamma, the angle between the MDC and Mdev on the 60-deg deviatoric arc
lammag = sqrt(lamdev(1,:).^2 + lamdev(2,:).^2 + lamdev(3,:).^2);
switch isort
    case 1
        gdot = (lamdev(1,:)-lamdev(3,:)) ./ (sqrt(2)*lammag);
        sg = sign(lamdev(2,:));
    case 2
        gdot = (lamdev(3,:)-lamdev(1,:)) ./ (sqrt(2)*lammag);
        sg = sign(lamdev(2,:));
    case 3
        gdot = (abs(lamdev(1,:))+abs(lamdev(2,:))) ./ (sqrt(2)*lammag);
        sg = sign(lamdev(3,:));
    case 4
        gdot = (abs(lamdev(3,:))+abs(lamdev(2,:))) ./ (sqrt(2)*lammag);
        sg = sign(lamdev(1,:));
end
gdot(gdot > 1) = 1; gdot(gdot <-1) = -1;
gamma = sg .* acos(gdot) * 180/pi;

% fraction of moment that is isotropic
fiso = zeros(n,1);
for ii=1:n
    % an invariant proxy for M0, which is never zero
    M0 = norm( lamdev(:,ii) + lamiso(:,ii) );
    fiso(ii) = norm(lamiso(:,ii)) / M0;
    fdev(ii) = norm(lamdev(:,ii)) / M0;
end
eISO = 100*fiso;
eDEV = 100*fdev;

% express as percent CLVD and percent DC -- Jost and Herrmann (1989)
eCLVD = 200*epsilon;
eDC   = 100*(1 - 2*abs(epsilon));

% If the tensor is purely isotropic, then there are no deviatoric
% eigenvalues, and we set epsilon, eCLVD, eDC all to zero.
iiso = find( abs(eISO-100) <= 1e-6 );
if ~isempty(iiso)
    eCLVD(iiso) = 0; eDC(iiso) = 0; epsilon(iiso) = 0;
end

% option to display figure (requires plot_histo.m)
ifigure = 0;
if ifigure==1
    nbin = 31; ymax = 0.12;
    if iepsilon==3, emax = 1.0; else emax = 0.5; end
    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot_histo(epsilon,linspace(-emax,emax,nbin)); ylim([0 ymax]);
    xlabel(sprintf('EPSILON (iepsilon = %i, isort = %i)',iepsilon,isort));
    subplot(nr,nc,2); plot_histo(gamma,linspace(-30,30,nbin)); ylim([0 ymax]);
    xlabel(sprintf('GAMMA, degrees (isort = %i)',isort));
end

%-------------------

if 0==1
    clear, clc, close all
    % get some test moment tensors from CMT catalog
    [otime,tshift,hdur,lat,lon,dep,M,~,~,eid,elabel,...
        str1,dip1,rk1,str2,dip2,rk2,lams,icmt,dataused,itri] = readCMT_all;
    n = length(otime);
    
    [gamma, epsilon, eCLVD, eDC, trM, eISO, eDEV] = CMT2epsilon(M,1,1);
    
    % check that changing isort will not affect the epsilon values (for fixed iepsilon)
    iepsilon = 1;
    epsall = zeros(n,4);
    gamall = zeros(n,4);
    for ii=1:4
        isort = ii;
        [gamma, epsilon, eCLVD, eDC, trM, eISO, eDEV] = CMT2epsilon(M,iepsilon,isort);
        epsall(:,ii) = epsilon;
        gamall(:,ii) = gamma;
    end
    figure; plot(epsall,'.');
    figure; plot(gamall,'.');

%     dir1 = '/net/sierra/raid1/carltape/results/SOURCES/socal_4/CMT_files_post_inverted/';
%     filename = [dir1 'CMTSOLUTION_9818433'];
%     [date,tshift,hdur,lat,lon,dep,M,eid,elabel] = readCMT(filename,13,0);
%     [epsilon, eCLVD, eDC, trM] = CMT2epsilon(M)
end

%==========================================================================
