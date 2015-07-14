function lam = nualpha2lam(nu,alpha)
%NUALPHA2LAM converts nu and alpha to unit lambda vector
%
% INPUT
%   nu      n-dimensional vector of Poisson ratios (unitless)
%   alpha   n-dimensional vector of angles between fault normal and slip vector, degrees [0,180]
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

% unit norm
lam = [lam1 ; lam2 ; lam3];
mag = sqrt( lam1.^2 + lam2.^2 + lam3.^2 );
lam = lam ./ repmat(mag,3,1);
%mag = sqrt( lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2 );  % check

%==========================================================================
% EXAMPLE

if 0==1
    %% grid of nu-alpha values
    avec = [5:10:175];           % avoid alpha = 90 (nu nonunique)
    nvec = [-1:0.05:0.45];       % restrict nu to permissible range
    [nu,alpha] = meshgrid(nvec,avec);
    nu = nu(:); alpha = alpha(:);
    lam = nualpha2lam(nu,alpha);
    
    % check
    [nucheck,alphacheck] = lam2nualpha(lam);
    norm(nucheck - nu)
    norm(alphacheck - alpha)
    
    %% CDC arc for nu=0.25
    avec = [0:2:180];
    lam = nualpha2lam(0.25,avec);
    [gamma,delta,M0,mu,lamdev,lamiso] = lam2lune(lam);
    figure; plot(gamma,delta,'.'); axis([-30 30 -90 90]);
end

%==========================================================================
