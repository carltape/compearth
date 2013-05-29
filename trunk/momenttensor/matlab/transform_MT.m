%
% function Mout = transform_MT(T,Min)
% Carl Tape, 11-March-2011
%
% This program will transform a set of symmetric matrices (e.g., moment
% tensors) using the transformation matrix T.
%
% moment tensor M = [Mrr Mtt Mpp Mrt Mrp Mtp]
%
% INPUT
%   T       3 x 3 transformation matrix
%   Min     6 x n input matrices: M = [M11 M22 M33 M12 M13 M23]
%
% OUTPUT
%   Mout    3 x 3 x n output matrices: M = [M11 M22 M33 M12 M13 M23]
% 
% calls transform_mat.m
% called by xxx
%

function Mout = transform_MT(T,Min)

% check that M is 6 x n
[Min,n] = Mdim(Min);

% transform to 3 x 3 x n
Min = Mvec2Mmat(Min,1);

% apply transformation T
Mout = transform_mat(T,Min);

% transform to 6 x n
Mout = Mvec2Mmat(Mout,0);

%=====================================================