%
% function [Aunit, Alen] = unit(A)
%
% This function inputs a 3 x n matrix of (x,y,z) points and ouputs the unit
% vectors, as well as the norms.
%
% EXAMPLE:
%   A = rand(3,6); [Aunit, Alen] = unit(A); repmat(Alen,3,1).*Aunit - A
%
% calls xxx
% called by xxx
%

function [Aunit, Alen] = unit(A)

% ensure that A is 3 x n
[m,n] = size(A);
if m ~= 3, A = A'; end

% Aunit = zeros(3,n);
% for i=1:n
%     Aunit(:,i) = A(:,i) / norm(A(:,i));
% end

% lengths of the n vectors (1 x n)
Alen = sqrt( A(1,:).^2 + A(2,:).^2 + A(3,:).^2 );

% unit vectors (3 x n)
Aunit = A ./ repmat(Alen,3,1);

%===================================
