function zhatn = zhatngammaomega(phi,sigma,gamma,omega)
% 
%

% note: no need to return b or c
[Z,W,a] = ZWgamma(phi,sigma,gamma);

% Eq 20a
zn = Z - sqrt( (cos(omega) - W)./a );
z = zn;

% Eq 20b
inds = find(cos(omega) < W);
z(inds) = Z(inds);

% Eq 26
zhatn = z;
zhatn(z < -1) = -1;
zhatn(z > 1) = 1;
