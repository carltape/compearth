function v = gamma2v(gamma)
%GAMMA2V v(gamma) to get v from lune longitude
%
% INPUT
%   gamma   n x 1 vector of gamma angles, radians (gamma = [-pi/6, pi/6])
%
% OUTPUT
%   v       n x 1 vector (v = [-1/3, 1/3])
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%
% Example values:
%   gamma2v(-pi/6) = -1/3
%   gamma2v(0) = 0
%   gamma2v(pi/6) = 1/3
%
% note: might want to check the input range of v
%

v = (1/3)*sin(3*gamma);
