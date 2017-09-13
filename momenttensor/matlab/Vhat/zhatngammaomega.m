function zhatn = zhatngammaomega(phi,sigma,gamma,omega)
%ZHATNGAMMAOMEGA functions needed for calculating Vhat_gamma(omega) curves
%
% W. Tape and C. Tape, 2017, GJI
% Volume in moment tensor space in terms of distance
%
% called by Vgammaomega.m
%

% note: no need to return b or c
[Z,W,a] = ZWgamma(phi,sigma,gamma);

% Tape and Tape (2017), Eq 20a
zn = Z - sqrt( (cos(omega) - W)./a );
z = zn;

% Tape and Tape (2017), Eq 20b
inds = find(cos(omega) < W);
z(inds) = Z(inds);

% Tape and Tape (2017), Eq 26
zhatn = z;
zhatn(z < -1) = -1;
zhatn(z > 1) = 1;
