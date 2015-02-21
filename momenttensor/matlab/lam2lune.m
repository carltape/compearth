function [gamma,delta,M0,thetadc,lamdev,lamiso] = lam2lune(lam)
%LAM2LUNE convert eigenvalues to lune coordinates (gamma, delta, M0)
%
% INPUT
%   lam         3 x n set of eigenvalues for a set of moment tensors
%               
% OUTPUT
%   gamma       angle from DC meridian to lune point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to lune point (-90 <= delta <= 90)
%   M0          seismic moment, M0 = ||lam|| / sqrt(2)
%   thetadc     angle from DC to lune point (0 <= thetadc <= 90)
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

[lam,n] = lamsort(lam);

% magnitude of lambda vector (rho of TT2012 -- see p. 490 within text)
%lammag = sqrt(2) * M0;
lammag = sqrt(lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2);

% seismic moment
%M0 = sqrt(lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2);
M0 = lammag / sqrt(2);

% TapeTape2012a, Eq. 21a (and 23)
% numerical safety 1: if trace(M) = 0, delta = 0
% numerical safety 2: is abs(bdot) > 1, adjust bdot to +1 or -1
delta = zeros(1,n);         % initialized to delta=0
idev = find(sum(lam) ~= 0);
bdot = (lam(1,:) + lam(2,:) + lam(3,:)) ./ (sqrt(3)*lammag);
bdot(bdot > 1) = 1; bdot(bdot <-1) = -1;
delta(idev) = 90 - acos(bdot(idev)) * deg;

% TapeTape2012a, Eq. 21a
% note: we set gamma=0 for (1,1,1) and (-1,-1,-1)
gamma = atan((-lam(1,:) + 2*lam(2,:) - lam(3,:)) ./ (sqrt(3)*(lam(1,:) - lam(3,:)))) * deg;
biso = lam(1,:)==lam(3,:);
gamma(biso) = 0;

% extra output
trM = sum(lam); 
lamiso = repmat(1/3*trM,3,1);
lamdev = lam - lamiso;

% compute thetadc -- the angle between the DC and the lune point
%theta = acos( cos(delta/deg) .* cos(gamma/deg) ) * deg;
thetadc = acos( (lam(1,:) - lam(3,:)) ./ (sqrt(2)*lammag) ) * deg;

% column vectors
delta = delta(:);
gamma = gamma(:);
M0 = M0(:);
thetadc = thetadc(:);

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
    lam = lune2lam(gamma0,delta0,M00);
    
    [gamma,delta,M0,thetadc] = lam2lune(lam);
    
    figure; nr=2; nc=2;
    subplot(nr,nc,1); plot(gamma-gamma0,'.'); title('gamma residual');
    subplot(nr,nc,2); plot(delta-delta0,'.'); title('delta residual');
    subplot(nr,nc,3); plot(M0-M00,'.'); title('M0 residual');
    
    figure; scatter(gamma,delta,4^2,thetadc,'filled');
    caxis([0 90]); colorbar;
    xlabel('gamma, deg'); ylabel('delta, deg');
    title('angle between DC and MT point');
end

%==========================================================================

