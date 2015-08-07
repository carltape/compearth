function v = gamma2v(gamma)
%GAMMA2V v(gamma) for lune longitude
%
% INPUT
%   gamma   n x 1 vector of gamma angles, radians
%
% OUTPUT
%   v       n x 1 vector
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%
% Example values:
%   gamma2v(-pi/6) = -1/3
%   gamma2v(0) = 0
%   gamma2v(pi/6) = 1/3
%

v = (1/3)*sin(3*gamma);
