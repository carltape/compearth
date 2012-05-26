function U = rotmat_gen(u,alpha)
% This returns a rotation matrix, given an angle and an index for the axis.
%
% INPUT
%   u       rotation axis
%   alpha   rotation angle, degrees
% OUTPUT
%   U       rotation matrix
%
% calls rotmat.m
%
% Carl Tape 5/2012
%

deg = 180/pi;

% get (phi,theta) for rotation axis
[uph,ele,rho] = cart2sph(u(1),u(2),u(3));
uth = pi/2 - ele;

% note: operations from right to left
U = rotmat(uph*deg,3)*rotmat(uth*deg,2)*rotmat(alpha,3)*rotmat(-uth*deg,2)*rotmat(-uph*deg,3);

%==========================================================================
