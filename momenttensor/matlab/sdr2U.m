function U = sdr2U(kappa,theta,sigma)
%SDR2U strike-dip-rake to rotation matrix U 
%
% INPUT
%   kappa   strike angle, degrees
%   theta   dip angle, degrees
%   sigma   slip (or rake) angle, degrees
%
% OUTPUT
%   U       3 x 3 x n set of bases in SOUTH-EAST-UP convention
%
% See WTape and CTape (2012) "A geometric setting for moment tensors" (TT2012).
%
% See examples in U2sdr.m
% called by TT2CMT.m
%
% Carl Tape, 2012/12
%

nk = length(kappa);
nt = length(theta);
ns = length(sigma);
if length(unique([nk nt ns])) > 2
    whos kappa theta sigma
    error('input must have same dimension');
end
n = nk;

for ii=1:n
    if theta(ii)==0
        warning('%i/%i input fault is horizontal, so strike angle (%.1f) is undefined',ii,n,kappa(ii));
        disp(sprintf('resetting slip angle (%.1f) to 0',sigma(ii)));
        sigma(ii) = 0;
    end
end

% NOTE: Algorithmically, it would be simpler to compute V directly from the
% expression in TT2012 Proposition 2, since this requires fewer calculations.
% (In the case here, we do not need the fault vectors at all.)
% The implementaion below is more conceptual.
% The basis is specified from the components of the north and zenith vectors.

% for north-west-up basis (TT2012)
%north = [1 0 0]'; zenith = [0 0 1]';

% for south-east-up basis (TT2013)
north = [-1 0 0]'; zenith = [0 0 1]';

% TT2012, p. 485
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
U = NaN(3,3,n);
Yrot = rotmat(-45,2);
for ii=1:n
    V = [S(:,ii) cross(N(:,ii),S(:,ii)) N(:,ii)];
    Ux = V*Yrot;
    U(:,:,ii) = Ux;
end

%==========================================================================
