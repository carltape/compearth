function mnorm = norm_MT(M,Lnorm)
%NORM_MT compute matrix (Frobenius) norm for symmetric matrix (moment tensor)
%
% INPUT
%   M       6 x n input symmetric matrices: M = [M11 M22 M33 M12 M13 M23]
%   Lnorm   optional: type of norm (default=2): 1,2,Inf,-Inf,etc
%
% OUTPUT
%   mnorm   n x 1 vector of matrix (Frobenius) norms
% 
% calls norm_mat.m
%
% Carl Tape, 11-March-2011
%

if nargin==1, Lnorm = 2; end

% check that M is 6 x n
[M,n] = Mdim(M);

% transform to 3 x 3 x n
M = Mvec2Mmat(M,1);

% compute norm
mnorm = norm_mat(M,Lnorm);

%==========================================================================