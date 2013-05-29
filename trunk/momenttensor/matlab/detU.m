function detU = detU(U)
%DETU compute determinant for a set of U matrices (U is 3 x 3 x n)

[a,b,n] = size(U);
if any([a b]~=3), error('U must be 3 x 3 x n'); end

detU = NaN(n,1);
for ii=1:n
   detU(ii) = det(squeeze(U(:,:,ii)));
end

%==========================================================================