function gamma = v2gamma(v)
%V2GAMMA gamma(v) to get lune longitude from v
%
% INPUT
%   v       n x 1 vector (v = [-1/3, 1/3])
%
% OUTPUT
%   gamma   n x 1 vector of gamma angles, radians (gamma = [-pi/6, pi/6])
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%
% Example values:
%   v2gamma(-1/3) = -pi/6
%   v2gamma(0) = 0
%   v2gamma(1/3) = pi/6
%
% note: might want to check the input range of v
%

gamma = (1/3)*asin(3*v);
