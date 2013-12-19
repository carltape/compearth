function [lam,n] = lamsort(lam)
%LAMSORT sort eigenvalues as lam1 >= lam2 >= lam3

[a,n] = size(lam);
if a~=3, error('dimension of lam (%i x %i) must be 3 x %i',a,n,n); end
if n==3, disp('warning: lam is 3 x 3, so make sure that each column corresponds to a lambda vector'); end

% our formulas assume that eigenvalues are sorted as lam1 >= lam2 >= lam3
lam = sort(lam,'descend');