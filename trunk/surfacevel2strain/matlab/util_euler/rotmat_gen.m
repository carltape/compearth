function U = rotmat_gen(v,xi)
%ROTMAT_GEN compute a rotation matrix, given an axis and an angle
%
% INPUT
%   v    rotation axis
%   xi   rotation angle, degrees
% OUTPUT
%   U    rotation matrix
%
% calls rotmat.m
%
% Carl Tape 5/2012

deg = 180/pi;

% get (phi,theta) for rotation axis
[vph,ele,rho] = cart2sph(v(1),v(2),v(3));
vth = pi/2 - ele;

% note: operations from right to left
U = rotmat(vph*deg,3)*rotmat(vth*deg,2)*rotmat(xi,3)*rotmat(-vth*deg,2)*rotmat(-vph*deg,3);

%==========================================================================
