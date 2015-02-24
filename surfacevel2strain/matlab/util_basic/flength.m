function m = flength(v)
%FLENGTH given a set of input vectors, compute their lengths
%
% INPUT
%   v      3 x n set of initial vectors
% OUTPUT
%   m      n x 1 set of vector lengths
%

[a,n] = size(v);

vdot = zeros(1,n);
for ii=1:a
   vdot = vdot + v(ii,:).^2;
end
m = sqrt( vdot );

m = m(:);

%----------------------

% % check dimensions
% [a,n] = size(v);
% if a ~= 3, error('v must be 3 x n'); end
% m = sqrt( v(1,:).^2 + v(2,:).^2 + v(3,:).^2 );
% m = m(:);
