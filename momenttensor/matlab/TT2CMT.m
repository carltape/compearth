function [M,lam,U] = TT2CMT(gamma,delta,M0,kappa,theta,sigma)
%TT2CMT convert geometrical parameters into moment tensors
%
% INPUT
%   gamma       angle from DC meridian to MT point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to MT point (-90 <= delta <= 90)
%                  note: delta = 0 is deviatoric
%   M0          seismic moment
%                  note: rho = sqrt(2)*M0, so set M0=1/sqrt(2) for rho=1
%   kappa       strike angle, degrees
%   theta       dip angle, degrees
%   sigma       slip (or rake) angle, degrees
%
% OUTPUT
%   M           6 x n set of moment tensors in CMT convention
%               M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%   lam         3 x n set of eigenvalues
%   U           3 x 3 x n set of bases in up-south-east convention
%
% Reverse program for CMT2TT.m
% See WTape and CTape (2012) "A geometric setting for moment tensors" (TT2012).
%
% Carl Tape, 12/2012
%

n = length(gamma);

for ii=1:n
    if theta(ii)==0
        disp(sprintf('WARNING: %i/%i input fault is horizontal, so strike angle (%.1f) is undefined',ii,n,kappa(ii)));
        disp(sprintf('         resetting slip angle (%.1f) to 0',sigma(ii)));
        sigma(ii) = 0;
    end
end

%---------------------
% PART 1: moment tensor source type (or pattern)

lam = lune2lam(gamma,delta,M0);

%---------------------
% PART 2: moment tensor orientation
% NOTE: Algorithmically, it would be simpler to compute V directly from the
% expression in Proposition 2. Here we do not need the fault vectors at
% all. But the implementaion below is more conceptual. If you want to
% return the fault vectors, then keep in mind the basis (north-west-up),
% or transform them to up-south-east to be consistent with M.

% convention in TT2012 is north-west-up
north = [1 0 0]';
zenith = [0 0 1]';

phi = -kappa;

% TT2012, Eq 27abc
K = NaN(3,n);
N = NaN(3,n);
S = NaN(3,n);
for ii=1:n
    K(:,ii) = rotmat(phi(ii),3) * north;
    N(:,ii) = rotmat_gen(K(:,ii),theta(ii)) * zenith;
    S(:,ii) = rotmat_gen(N(:,ii),sigma(ii)) * K(:,ii);
end

% TT2012, Eq 28 (or Proposition 2)
M = NaN(3,3,n);
U = NaN(3,3,n);
Yrot = rotmat(-45,2);
for ii=1:n
    V = [S(:,ii) cross(N(:,ii),S(:,ii)) N(:,ii)];
    Ux = V*Yrot;
    M(:,:,ii) = Ux*diag(lam(:,ii))*Ux';
    U(:,:,ii) = Ux;
end
M = Mvec2Mmat(M,2);

% convert from north-west-up to up-south-east
i1 = 3;
i2 = 1;
M = convert_MT(i1,i2,M);
U = convertv(i1,i2,U);
    
%==========================================================================
% EXAMPLES

if 0==1
    M0 = 1;  % seismic moment
    kappa = 320; theta = 10; sigma = 20;
    M = TT2CMT(0,0,M0,kappa,theta,sigma)
    [gamma,delta,M0,kappa,theta,sigma] = CMT2TT(M)
    
    % compare with Mathematica
    deg = 180/pi;
    theta = (pi/4)*deg;
    kappa = (pi/3)*deg;
    sigma = -(pi/3)*deg;
    rho = 1;
    gamma = (pi/12)*deg;
    beta = (pi/3)*deg;
    M0 = rho/sqrt(2);
    delta = 90-beta;
    M = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
    Maki = Mvec2Mmat(convert_MT(1,2,M),1)
end

%==========================================================================
