function [gamma,delta,M0,mu,lamdev,lamiso] = lam2lune(lam)
%LAM2LUNE convert eigenvalues to lune coordinates (gamma, delta, M0)
%
% INPUT
%   D           3 x n set of input moment tensors in eigenbasis
%               note 1: eigenvalues sorted as lam1 >= lam2 >= lam3
%               note 2: normalized such that each MT has moment M0
%               
% OUTPUT
%   gamma       angle from DC meridian to MT point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to MT point (-90 <= delta <= 90)
%   M0
%   mu
%   lamdev      eigenvalues of deviatoric component
%   lamiso      eigenvalues of isotropic component
%
% Reverse program for lune2lam.m
% See also CMT2all.m
%
% See TapeTape2012 "A geometric setting for moment tensors".
%
% Carl Tape, 01-April-2011
%

deg = 180/pi;

[a,n] = size(lam);
if a~=3, error('dimension of lam (%i x %i) must be 3 x %i',a,n,n); end

% magnitude of lambda vector (rho of TT2012 -- see p. 490 text)
%lammag = sqrt(2) * M0;
lammag = sqrt(lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2);

% seismic moment
%M0 = sqrt(lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2);
M0 = lammag / sqrt(2);

% decompose into isotropic and deviatoric parts
% note: this is not needed for the calculations below
trM = sum(lam); 
lamiso = repmat(1/3*trM,3,1);
lamdev = lam - lamiso;

% % gamma (use deviatoric eigenvalues)
% gamma = zeros(n,1);
% lamdevmag = sqrt(lamdev(1,:).^2 + lamdev(2,:).^2 + lamdev(3,:).^2);
% inoniso = find(lamdevmag > 0);
% if ~isempty(inoniso)
%     gdot(inoniso) = (lamdev(1,inoniso)-lamdev(3,inoniso)) ./ (sqrt(2)*lamdevmag(inoniso));
%     sg = sign(lamdev(2,:));
%     gdot(gdot > 1) = 1; gdot(gdot <-1) = -1;
%     gamma(inoniso) = sg(inoniso) .* acos(gdot(inoniso)) * deg;
% end

% TapeTape2012a, Eq. 21a (and 23)
% compute gamma and delta using ACTUAL eigenvalues (allowing for isotropic component)
% numerical safety 1: if trace(M) = 0, delta = 0
% numerical safety 2: is abs(bdot) > 1, adjust bdot to +1 or -1
delta = zeros(n,1);
idelta = find(sum(lam) ~= 0);
bdot = (lam(1,:) + lam(2,:) + lam(3,:)) ./ (sqrt(3)*lammag);
bdot(bdot > 1) = 1; bdot(bdot <-1) = -1;
delta(idelta) = 90 - acos(bdot(idelta)) * deg;

% TapeTape2012a, Eq. 21a
gamma = atan((-lam(1,:) + 2*lam(2,:) - lam(3,:)) ./ (sqrt(3)*(lam(1,:) - lam(3,:)))) * deg;

% column vectors
delta = delta(:);
gamma = gamma(:);
M0 = M0(:);

% conpute mu -- the angle between the DC (0,0) and the point
%mu = acos( cos(delta/deg) .* cos(gamma/deg) ) * deg;
mu = acos( (lam(1,:) - lam(3,:)) ./ (sqrt(2)*lammag) ) * deg;

%==========================================================================
% EXAMPLE

if 0==1
    clear, close all, clc
    gvec = linspace(-30,30,100);
    bvec = linspace(-89,89,100);
    [G,B] = meshgrid(gvec,bvec);
    gamma0 = G(:);
    delta0 = B(:);
    M00 = 1e16*ones(length(gamma0),1);
    D = lune2lam(gamma0,delta0,M00);
    
    [gamma,delta,M0,mu] = lam2lune(D);
    
    figure; nr=2; nc=2;
    subplot(nr,nc,1); plot(gamma-gamma0,'.'); title('gamma residual');
    subplot(nr,nc,2); plot(delta-delta0,'.'); title('delta residual');
    subplot(nr,nc,3); plot(M0-M00,'.'); title('M0 residual');
    
    figure; scatter(gamma,delta,4^2,mu,'filled');
    caxis([0 90]); colorbar;
    xlabel('gamma, deg'); ylabel('delta, deg');
    title('angle between DC and MT point');
end

%==========================================================================

