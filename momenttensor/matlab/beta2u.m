function u = beta2u(beta)
%BETA2V u(beta) for lune colatitude
%
% INPUT
%   beta    n x 1 vector of lune colatitudes, radians
%
% OUTPUT
%   u       n x 1 vector
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%

u = (3/4)*beta - (1/2)*sin(2*beta) + (1/16)*sin(4*beta);
