function shatn = shatngammaomega(mu,nu,gamma,omega)
% 
%

[a,b,c] = abcgamma(mu,nu,gamma);
S = -b ./ (2*a);
T = -(b.^2 - 4*a.*c)./(4*a);
temp = sqrt(b.^2 - 4*a.*(c - cos(omega)));
rm = (-b - temp) ./ (2*a);

inds = find(cos(omega) < T);
s = rm;
s(inds) = S(inds);

shatn = s;
shatn(s < 0) = 0;
shatn(s > 1) = 1;
