function U = UiU(U1,U2)
%UiU input U1 and U2 and output U1^-1 * U2
%
% INPUT
%   U1     3 x 3 x n set of rotation matrices
%   U2     3 x 3 x n set of rotation matrices
%
% OUTPUT   
%   UiU    n x 1 set of rotation angles
%

% check dimensions
[~,~,n1] = size(U1);
[~,~,n2] = size(U2);
if n1~=n2, error('U1 and U2 must be equal in dimension'); end
n = n1;

% compute U = U1' * U2
U = NaN(size(U1));
for ii=1:n
    U(:,:,ii) = U1(:,:,ii)' * U2(:,:,ii);
end
