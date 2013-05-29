function [f1,f2,f3] = fmat2fvec(fmat)
% extract fault vectors from a matrix of values

[a,n] = size(fmat);
if a~=9, error('fmat must be 9 x n'); end

f1 = fmat(1:3,:);
f2 = fmat(4:6,:);
f3 = fmat(7:9,:);
