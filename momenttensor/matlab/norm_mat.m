function mnorm = norm_mat(Min,Lnorm)
%NORM_MAT compute matrix (Frobenius) norms for a set of input matrices
%
% INPUT
%   Min     3 x 3 x n input matrices
%   Lnorm   optional: type of norm (default=2): 1,2,Inf,-Inf,etc
%
% OUTPUT
%   mnorm   n x 1 vector of matrix (Frobenius) norms
% 
% called by norm_MT.m
%
% Carl Tape 11/2010

if nargin==1, Lnorm = 2; end

% check that M is 3 x 3 x n
[a,b,n] = size(Min);
if any([a b]~=3), error('M should be 3 x 3 x n'); end

% resize to 9 x n
M = reshape(Min,9,n);

mnorm = zeros(n,1);
for ii = 1:n
    mnorm(ii) = norm(M(:,ii),Lnorm);
end

%==========================================================================