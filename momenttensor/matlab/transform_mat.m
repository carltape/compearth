function Mout = transform_mat(T,Min)
%TRANSFORM_MT transform a set of 3 x 3 tensors using transformation matrix T
%
% INPUT
%   T       3 x 3 transformation matrix
%   Min     3 x 3 x n input matrices
%
% OUTPUT
%   Mout    3 x 3 x n output matrices
%
% The version associated with symmetric M is transform_MT.m
%
% Carl Tape, 17-Mar-2011
%

% check that M is 3 x 3 x n
[a,b,n] = size(Min);
if any([a b]~=3), error('M should be 3 x 3 x n'); end

Mout = 0*Min;
for ii = 1:n
    M0 = Min(:,:,ii);
    Mout(:,:,ii) = T * M0 * T';
end
