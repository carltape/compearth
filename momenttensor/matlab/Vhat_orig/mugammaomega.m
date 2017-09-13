function mugamma = mugammaomega(gamma,omega)
% 
%
% INPUT
%   gamma   in radians
%   omega   in radians
%
% called by Vgammaomega.m
%

% calculate mu_gamma(omega)
n = length(gamma);
gtemp = pi/3 - 2*gamma;
inds = find( omega <= gtemp );
mugamma = pi/2 * ones(n,1);
mugamma(inds) = asin( sin(omega(inds)/2) ./ sin(pi/6 - gamma(inds)) );
