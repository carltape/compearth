function fmat = fvec2fmat(f1,f2,f3)
% conver from a set of fault vectors to a matrix of numbers

[a1,n1] = size(f1);
[a2,n2] = size(f2);
[a3,n3] = size(f3);
if ~and( length(unique([n1 n2 n3])==1), all([a1 a2 a3]==3))
    error('f1, f2, f3 must have the same lengths');
end

fmat = [f1 ; f2 ; f3];
