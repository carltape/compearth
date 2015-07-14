function [phi,zeta] = lam2phizeta(lam)
%LAM2PHIZETA converts eigenvalues to phi and zeta
%
% Convert eigenvalues of a moment tensor into two quantities that represent
% the crack-plus-double-couple (CDC) model, phi and zeta.
%
% INPUT
%   lam     3 x n set of eigenvalue triples
%
% OUTPUT
%   phi     n x 1 vector of phi angles on the lune, degrees [-180,180]
%   zeta    n x 1 vector of crack fraction within CDC model [0,90]
%
% Reverse function is phizeta2lam.m
% See Tape and Tape (2013), "The classical model for moment tensors"
% 
% Carl Tape, 08-Jan-2013
%

deg = 180/pi;

lam = lamsort(lam);

lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);

rho = sqrt(sum(lam.^2));

% TT2013, Eqs 24
phi = atan2( lam1-2*lam2+lam3, sqrt(2)*(lam1+lam2+lam3) )*deg;

zeta = acos( sqrt(2*(lam1-lam2).*(lam2-lam3)) ./ rho )*deg;

% column vectors
phi = phi(:);
zeta = zeta(:);

%==========================================================================
% EXAMPLES

if 0==1
    %% example from TT2013, App A.
    lam = [ 8.802 2.584 -1.851]';
    [phi,zeta] = lam2phizeta(lam)
    lamcheck = phizeta2lam(phi,zeta)
    lam / norm(lam)
    
    %% grid of values on the lune
    dd = 0.5;
    gvec = [-29:dd:29];
    dvec = [-89:dd:89];
    [Gamma,Delta] = meshgrid(gvec,dvec);
    gamma = Gamma(:); delta = Delta(:);
    lam = lune2lam(gamma,delta);
    % calculate phi and zeta
    [phi,zeta] = lam2phizeta(lam);
    % scatter plot
    figure; subplot(1,2,1); scatter(gamma,delta,4^2,phi,'filled');
    axis equal, axis tight; caxis([-180 180]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\phi');
    subplot(1,2,2); scatter(gamma,delta,4^2,zeta,'filled');
    axis equal, axis tight; caxis([0 90]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\zeta');
    % contour lines
    phicons = [-180:10:180];
    figure; subplot(1,2,1); contour(Gamma,Delta,reshape(phi,size(Gamma)),phicons);
    axis equal, axis tight; caxis([-180 180]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\phi');
    zcons = [0:10:90];
    subplot(1,2,2); contour(Gamma,Delta,reshape(zeta,size(Gamma)),zcons);
    axis equal, axis tight; caxis([0 90]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\zeta');
end
    
%==========================================================================
