function [Aunit, Alen] = unit(A)
%UNIT output unit vectors, given a 3 x n input matrix of vectors
%
% EXAMPLE:
%   A = rand(3,6); [Aunit, Alen] = unit(A); repmat(Alen,3,1).*Aunit - A
%
% JUST USE normc OR normr INSTEAD OF THIS FUNCTION
%

% ensure that A is 3 x n
[m,~] = size(A);
if m ~= 3, A = A'; end

% lengths of the n vectors (1 x n)
Alen = sqrt( A(1,:).^2 + A(2,:).^2 + A(3,:).^2 );

% unit vectors (3 x n)
Aunit = A ./ repmat(Alen,3,1);
