%
% function R = euler2rotmat(evec)
% Carl Tape, 01-Nov-2005
%
% Given an euler pole with corresponding angular rotation angle, this
% program returns the rotation matrix that is to be applied to a set of
% Cartesian (x,y,z) points on a sphere.
%
% This also allows for "non-standard" latitude, longitude, and rotation
% angles, and returns the "standard" form.  For example, an euler vector
% specified by (latitude = -120, longitude = 465, and omega = -300) is
% permissible, but the "standard" version will be returned
%
% INPUT: euler pole vector
%    evec(1) = latitude (deg) of euler pole
%    evec(2) = longitude (deg) of euler pole
%    evec(3) = rotation angle (deg)
%
% OUTPUT:
%    R       = matrix for finite rotation
%    evec    = euler pole vector (lat, lon=[-180,180], omeg > 0)
%
% Formulas from Cox and Hart, Plate Tectonics (1986), p. 226.
% EXAMPLE: R = euler2rotmat([-37 312 65])
%
% Reverse program is rotmat2euler.m
%
% See eulerREADME for related programs.
%
% calls euler_mapping.m, latlon2xyz.m
% called by euler_rot_tec.m, test_euler_rot_tec.m
%

function [R,evec] = euler2rotmat(evec)

deg = 180/pi;

% euler vector
elat = evec(1);
elon = evec(2);
omega = evec(3);

% shift to "standard coordinates", including a positive rotation angle
evec = euler_mapping(evec);

% unit vector for rotation pole
exyz = latlon2xyz(elat,elon,1);
ex = exyz(1);
ey = exyz(2);
ez = exyz(3);

% rotation matrix
R = zeros(3,3);
coso = cos(omega/deg);
sino = sin(omega/deg);
R(1,1) = ex*ex*(1 - coso) + coso;
R(1,2) = ex*ey*(1 - coso) - ez*sino;
R(1,3) = ex*ez*(1 - coso) + ey*sino;
R(2,1) = ey*ex*(1 - coso) + ez*sino;
R(2,2) = ey*ey*(1 - coso) + coso;
R(2,3) = ey*ez*(1 - coso) - ex*sino;
R(3,1) = ez*ex*(1 - coso) - ey*sino;
R(3,2) = ez*ey*(1 - coso) + ex*sino;
R(3,3) = ez*ez*(1 - coso) + coso;

%==============================================================
