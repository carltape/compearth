function gamma = v2gamma(v)
%V2GAMMA v(gamma) for lune longitude
%
% OUTPUT
%   gamma   n x 1 vector of gamma angles, radians
%
% INPUT
%   v       n x 1 vector
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%

gamma = (1/3)*asin(3*v);
