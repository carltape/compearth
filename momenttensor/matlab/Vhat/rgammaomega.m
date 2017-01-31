function [rp,rm,a,b,c] = rgammaomega(mu,nu,gamma,omega)
% 
%
% INPUT
%   mu      in radians
%   nu      in radians
%   gamma   in radians
%   omega   in radians
%

[a,b,c] = abcgamma(mu,nu,gamma);

temp = sqrt(b.^2 - 4*a.*(c - cos(omega)));
rp = (-b + temp) ./ (2*a);
rm = (-b - temp) ./ (2*a);
