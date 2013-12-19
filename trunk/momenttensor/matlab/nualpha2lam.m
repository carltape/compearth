function lam = nualpha2lam(nu,alpha)
%NUALPHA2LAM converts nu and alpha to unit lambda vector
%
% INPUT
%   nu      n x 1 vector of Poisson ratios (unitless)
%   alpha   n x 1 vector of angles between fault normal and slip vector, degrees
%
%  OUTPUT
%   lam     3 x n set of normalized eigenvalue triples, sorted lam1 >= lam2 >= lam3
%
% Reverse function is lam2nualpha.m
% See Tape and Tape (2013), "The classical model for moment tensors"
% 
% Carl Tape, 08-Jan-2013
%

% row vectors
alpha = alpha(:)';
nu = nu(:)';

cosa = cos(alpha*pi/180);

% TapeTape2013, Eq 30
lam1 = cosa ./ (1-2*nu) + 1;
lam2 = 2 .* nu .* cosa ./ (1-2*nu);
lam3 = cosa ./ (1-2*nu) - 1;

lam = [lam1 ; lam2 ; lam3];
mag = sqrt( lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2 );
lam = lam ./ repmat(mag,3,1);
%mag = sqrt( lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2 );  % check

%==========================================================================
% EXAMPLE

if 0==1
    % grid of nu-alpha values
    numa = 10;
    numn = numa;
    avec = linspace(0,180,numa);
    nvec = linspace(-1,0.49,numn);
    [alpha,nu] = meshgrid(avec,nvec);
    alpha = alpha(:);
    nu = nu(:);
    lam = nualpha2lam(nu,alpha);
    
    % CDC arc for nu=0.25
    numa = 100;
    avec = linspace(0,180,numa);
    lam = nualpha2lam(0.25,avec);
    [gamma,delta,M0,mu,lamdev,lamiso] = lam2lune(lam);
    figure; plot(gamma,delta,'.'); axis([-30 30 -90 90]);
end

%==========================================================================
