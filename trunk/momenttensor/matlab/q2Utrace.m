function trU = q2Utrace(q)
%Q2UTRACE convert unit quaternion to trace(U)
%
% INPUT   
%   q       4 x n set of unit quaternions
%
% OUTPUT
%   trU     1 x n set of tr(U)
% 
% See q2U.m, which will also return U.
%
% Carl Tape, 9/26/2013
%

w = q(1,:);
x = q(2,:);
y = q(3,:);
z = q(4,:);
trU = 3*w.^2 - x.^2 - y.^2 - z.^2;
