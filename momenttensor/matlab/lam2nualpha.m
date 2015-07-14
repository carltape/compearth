function [nu,alpha] = lam2nualpha(lam)
%LAM2NUALPHA converts eigenvalues to nu and alpha
%
% INPUT
%   lam     3 x n set of eigenvalue triples
%
% OUTPUT
%   nu      n x 1 vector of Poisson parameters (unitless)
%   alpha   n x 1 vector of angles between fault normal and slip vector, degrees [0,180]
%
% Reverse function is nualpha2lam.m
% See Tape and Tape (2013), "The classical model for moment tensors"
% 
% Note that nu can take on any value for the double couple (alpha = 90).
%
% Carl Tape, 08-Jan-2013
%

lam = lamsort(lam);

lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);

% TT2013, Eqs 32ab
alpha = 180/pi* acos( (lam1 - 2*lam2 + lam3) ./ (lam1 - lam3) );
nu = lam2 ./ (lam1 + lam3);

% column vectors
alpha = alpha(:);
nu = nu(:);

%==========================================================================
% EXAMPLES

if 0==1
    %% example from TT2013, App A.
    lam = [ 8.802 2.584 -1.851]';
    [nu,alpha] = lam2nualpha(lam)
    lamcheck = nualpha2lam(nu,alpha)
    lam / norm(lam)
    
    %% grid of values on the lune
    dd = 0.5;
    gvec = [-29:dd:29];
    dvec = [-89:dd:89];
    [Gamma,Delta] = meshgrid(gvec,dvec);
    gamma = Gamma(:); delta = Delta(:);
    lam = lune2lam(gamma,delta);
    % calculate nu and alpha
    [nu,alpha] = lam2nualpha(lam);
    % scatter plot
    figure; subplot(1,2,1); scatter(gamma,delta,4^2,nu,'filled');
    axis equal, axis tight; caxis([-1 0.5]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\nu');
    subplot(1,2,2); scatter(gamma,delta,4^2,alpha,'filled');
    axis equal, axis tight; caxis([0 180]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\alpha');
    % contour lines
    nucons = [-2 -1:0.1:0.5 2];
    figure; subplot(1,2,1); contour(Gamma,Delta,reshape(nu,size(Gamma)),nucons);
    axis equal, axis tight; caxis([-1 0.5]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\nu');
    subplot(1,2,2); contour(Gamma,Delta,reshape(alpha,size(Gamma)));
    axis equal, axis tight; caxis([0 180]); colorbar;
    xlabel('lune longitude'); ylabel('lune latitude'); title('\alpha');
end
    
%==========================================================================
