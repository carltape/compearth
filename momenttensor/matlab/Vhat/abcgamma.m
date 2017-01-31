function [a,b,c] = abcgamma(mu,nu,gamma)
% 
%
% INPUT
%   mu      in radians
%   nu      in radians
%   gamma   in radians
%

cosmu  = cos(mu);
cosnu  = cos(nu);
cos2mu = cos(2*mu);
cos2nu = cos(2*nu);
cos2gam = cos(2*gamma);
sin2gam = sin(2*gamma);

a = 1/16       * (6 - cos2mu - 12*cosmu.*cosnu - cos2nu).*cos2gam ...
    +sqrt(3)/8 * (4 - (cosmu - cosnu).^2).*sin2gam + ...
    +1/8       * (6 + cos2mu + cos2nu);

b = 1/2 * (cos2mu - cos2nu) .* sin(pi/6 - gamma).^2;

c = -1/16      * (2 + cos2mu - 12*cosmu.*cosnu + cos2nu) .* cos2gam ...
    -sqrt(3)/8 * (cosmu + cosnu).^2 .* sin2gam ...
    -1/8       * (2 - cos2mu - cos2nu);
