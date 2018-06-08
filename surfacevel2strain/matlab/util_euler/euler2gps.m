function [Vrtp, Vxyz] = euler2gps(evec, Pxyz)
%EULER2GPS calculate surface velocities from input points and euler vectors
%
% This function computes the linear velocity (vx,vy,vz) due to angular rotation
% at a set of input points P(x,y,z), which are expressed in METERS.
% See Cox and Hart (1986), p. 155.
%
% V = (w E) x (r P), where E is the Euler pole unit vector, expressed as (wx,wy,wz),
% where omega = sqrt(wx^2 + wy^2 + wz^2) and is in UNITS OF DEG/MYR.
%
% To convert to local coordinates, use global2local.m (see plate_model.m).
%
% INPUT:
%    evec   euler pole (wx,wy,wz), omega = norm(evec)
%    Pxyz   points in space (meters), which revolve around evec
% OUTPUT:
%    Vrtp   velocities in local basis (mm/yr)
%    Vxyz   velocities in global basis (mm/yr)
%
% Reverse program is gps2euler.m
%
% calls global2local.m, unit.m
% called by plate_model.m, test_gps2euler.m
%
% Carl Tape, 2005-10-21
%

% euler pole : convert deg/Myr --> rad/yr
% NOTE: this is not a unit vector
evec = evec * 1e-6 * pi/180;
ex = evec(1);
ey = evec(2);
ez = evec(3);

% ensure that Pxyz is 3 x n
[n,m] = size(Pxyz); if n~=3, Pxyz = Pxyz'; end

% local observation points : convert m --> mm
% EXAMPLE : norm(Pxyz(:,ii)) = 6.3710e+09 mm
% NOTE: this is not a unit vector
Pxyz = Pxyz * 1e3;
px = Pxyz(1,:);
py = Pxyz(2,:);
pz = Pxyz(3,:);

% V = (w E) x (r P), where E is the Euler pole unit vector
% if omega == 0, then the surace velocities =0
Vxyz = zeros(3,length(px));
Vrtp = zeros(3,length(px));
if norm(evec) ~= 0
    Vxyz(1,:) = ey*pz - py*ez;
    Vxyz(2,:) = ez*px - pz*ex;
    Vxyz(3,:) = ex*py - px*ey;
    
    Vrtp = global2local(Vxyz, Pxyz, 1);
end

% disp('E'); evec
% disp('P'); Pxyz(:,1)
% disp('V'); Vxyz(:,1)
% disp('err'); mean(abs( cross( evec , Pxyz(:,1) ) - Vxyz(:,1) ))

%=============================

% % euler pole : convert deg/Myr --> rad/yr
% omega = norm(evec) * 1e-6 * pi/180;
% u_evec = evec / norm(evec);
% ex = u_evec(1);
% ey = u_evec(2);
% ez = u_evec(3);
% 
% % ensure that Pxyz is 3 x n
% [n,m] = size(Pxyz); if n~=3, Pxyz = Pxyz'; end
% 
% % get the lengths of P : convert m --> mm
% % (in general, these should all be on a sphere with uniform radius)
% [u_Pxyz, rvec] = unit(Pxyz);
% rvec = rvec(:)' * 1e3;         % row vector
% 
% % unit vector components
% px = u_Pxyz(1,:);
% py = u_Pxyz(2,:);
% pz = u_Pxyz(3,:);
% 
% % V = (w E) x (r P), where E is the Euler pole unit vector
% % if omega == 0, then the surace velocities =0
% Vxyz = zeros(3,length(px));
% if omega ~= 0
%     Vxyz(1,:) = omega * rvec .* (ey*pz - py*ez);
%     Vxyz(2,:) = omega * rvec .* (ez*px - pz*ex);
%     Vxyz(3,:) = omega * rvec .* (ex*py - px*ey);
% end
% 
% Vrtp = global2local(Vxyz, Pxyz, 1);

%==============================================================
