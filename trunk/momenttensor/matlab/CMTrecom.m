function M = CMTrecom(lam,U)
%CMTRECOM combine basis and eigenvalues to get moment tensor
%
% INPUT
%   lam     3 x n set of eigenvalues
%   U       3 x 3 x n set of bases
%
% OUTPUT
%   M       6 x n moment tensors in CMT convention
%           M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%
% Inverse program to CMTdecom.m
% See also transform_MT.m
%
% Carl Tape, 01-Nov-2007
%

%disp('CMTrecom.m:'); whos lam Q

[~,n1] = size(lam);
[~,~,n2] = size(U);
if n1 ~= n2
    n1, n2
    error('n1 and n2 must be the same');
else
    n = n1;
end

Marray = zeros(3,3,n);
for ii = 1:n
    U0 = U(:,:,ii);
    D0 = diag( lam(:,ii) );
    Marray(:,:,ii) = U0 * D0 * U0';
end

M = Mvec2Mmat(Marray,0);

%==========================================================================