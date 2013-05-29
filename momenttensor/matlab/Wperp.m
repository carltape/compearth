%
% function W = Wperp(U)
% Carl Tape, 12-April-2011
%
% INPUT
%   U       3 x n set of eigenvectors of rotation matrix associated with eigenvalue 1
% OUTPUT
%   W       3 x n set of rotation vectors
% 
% calls xxx
% called by xxx
%

function W = Wperp(U)

[a,n] = size(U);
if a~=3, error('U must be 3 x n'); end

W = zeros(3,n);
for ii=1:n
    u0 = U(:,ii);
    if all(cross(u0, [0 0 1]') == 0)
        W(:,ii) = cross(u0, [1 0 0]');
    else
        W(:,ii) = cross(u0, [0 0 1]');
    end
end

%=====================================================