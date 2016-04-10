function U = rotmat_gen(v,xi,bdisplay)
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
if 1==1
    [vph,ele,rho] = cart2sph(v(1),v(2),v(3));
    vth = pi/2 - ele;
else
    rho = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
    vth = acos(v(3) / rho);
    vph = atan2(v(2),v(1));
end

% note: operations from right to left
R1 = rotmat(-vph*deg,3);
R2 = rotmat(-vth*deg,2);
R3 = rotmat(xi,3);
R4 = rotmat(vth*deg,2);
R5 = rotmat(vph*deg,3);
U = R5*R4*R3*R2*R1;

if nargin==3
    if bdisplay
    disp('----------------');
    disp('rotmat_gen.m');
    disp(sprintf(' rotation axis is v = (%.3f, %.3f, %.3f)',v(1),v(2),v(3)));
    disp(sprintf('            v/||v|| = (%.3f, %.3f, %.3f)',v(1)/norm(v),v(2)/norm(v),v(3)/norm(v)));
    disp(sprintf('  polar angle theta = %.3f deg (%.3f rad)',vth*deg,vth));
    disp(sprintf('azimuthal angle phi = %.3f deg (%.3f rad)',vph*deg,vph));
    R1,R2,R3,R4,R5
    U
    disp('----------------');
    end
end

%==========================================================================
