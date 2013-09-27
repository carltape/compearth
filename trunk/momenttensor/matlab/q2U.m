function [U,trU] = q2U(q)
%Q2U convert unit quaternion to rotation matrix
%
% INPUT   
%   q       4 x n set of quaternions
%
% OUTPUT
%   U       3 x 3 x n set of rotation matrices
%   trU     1 x n set of tr(U)
%
% See q2Utrace.m if you only want the trace.
% 
% Carl Tape, 8/12/2012
%

[~,n] = size(q);
U = NaN(3,3,n);

% convert quaternions to exact unit quaternions
for ii=1:n
   q(:,ii) = q(:,ii) / norm(q(:,ii));
end

w = q(1,:);
x = q(2,:);
y = q(3,:);
z = q(4,:);

% TapeTape2012, Eq 7
U(1,1,:) = w.^2 + x.^2 - y.^2 -z.^2;
U(1,2,:) = 2*(x.*y - w.*z);
U(1,3,:) = 2*(w.*y + x.*z);
U(2,1,:) = 2*(x.*y + w.*z);
U(2,2,:) = w.^2 - x.^2 + y.^2 -z.^2;
U(2,3,:) = 2*(y.*z - w.*x);
U(3,1,:) = 2*(x.*z - w.*y);
U(3,2,:) = 2*(w.*x + y.*z);
U(3,3,:) = w.^2 - x.^2 - y.^2 + z.^2;

% trace of U
trU = 3*w.^2 - x.^2 - y.^2 - z.^2;

%==========================================================================
% EXAMPLES

if 0==1
    % testing orthogonalization
    clear, close all, clc, format long
    P = [-1 0 0 ; 0 -1 0 ; 0 0 1];
    iorthoU = 0;
    U1 = P * U2pa([24 120 41 232],0,iorthoU);
    disp('INPUT U1:');
    [xi0,xi,q,ixi0,trALL,imaxtr] = U2xi0(U1,2,1);
    U1
    disp('U1o3 (quaternions):');
    U1o3 = q2U(q)
    U1o3'*U1o3
end

%==========================================================================
