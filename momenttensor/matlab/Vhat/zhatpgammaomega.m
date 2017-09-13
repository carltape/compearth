function zhatp = zhatpgammaomega(phi,sigma,gamma,omega)
%ZHATPGAMMAOMEGA functions needed for calculating Vhat_gamma(omega) curves
%
% W. Tape and C. Tape, 2017, GJI
% Volume in moment tensor space in terms of distance
%
% called by Vgammaomega.m
%

% note: no need to return b or c
[Z,W,a] = ZWgamma(phi,sigma,gamma);

% Tape and Tape (2017), Eq 20a
zp = Z + sqrt( (cos(omega) - W)./a );
z = zp;

% Tape and Tape (2017), Eq 20b
inds = find(cos(omega) < W);
z(inds) = Z(inds);

% Tape and Tape (2017), Eq 26
zhatp = z;
zhatp(z < -1) = -1;
zhatp(z > 1) = 1;
