function zhatp = zhatpgammaomega(phi,sigma,gamma,omega)
% 
%

% note: no need to return b or c
[Z,W,a] = ZWgamma(phi,sigma,gamma);

% Eq 20a
zp = Z + sqrt( (cos(omega) - W)./a );
z = zp;

% Eq 20b
inds = find(cos(omega) < W);
z(inds) = Z(inds);

% Eq 26
zhatp = z;
zhatp(z < -1) = -1;
zhatp(z > 1) = 1;
