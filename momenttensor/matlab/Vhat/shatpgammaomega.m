function shatp = shatpgammaomega(mu,nu,gamma,omega)
% 
%

[a,b,c] = abcgamma(mu,nu,gamma);
S = -b ./ (2*a);
T = -(b.^2 - 4*a.*c)./(4*a);
temp = sqrt(b.^2 - 4*a.*(c - cos(omega)));
rp = (-b + temp) ./ (2*a);

inds = find(cos(omega) < T);
s = rp;
s(inds) = S(inds);

shatp = s;
shatp(s < 0) = 0;
shatp(s > 1) = 1;
