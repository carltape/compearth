%
% function [evec, ve_est, vn_est, Pxyz] = gps2euler(dlon, dlat, drad, ve, vn, se, sn)
% Pablo Muse and Carl Tape, 21-July-2008
%
% This is a non-trivial program, unlike the reverse program euler2gps.m.
% For V = E x P, given V and P, there is no unique E.
% Hence we need a least squares solution to determine E, assuming all the V
% were derived from the same Euler pole (i.e., rotation vector).
%
% INPUT:
%    ve, vn: Local velocity field, M/YR
%    se, sn: Uncertainties in velocity components, M/YR
%     ---> THIS IS NOT PRESENTLY IMPLEMENTED
%
% OUTPUT:
%    evec: Best-fitting euler pole, [elat elon omega] with omega in deg/Myr
%    Pxyz: observations points in Cartesian coordinates (xyz), in METERS
%    ve_est, vn_est: estimated local velocities, in M/YR
%
% EXAMPLE:
%
% calls global2local.m, rotmat2euler.m
% called by test_gps2euler.m
%

function [evec, ve_est, vn_est, Pxyz] = gps2euler(dlon, dlat, drad, ve, vn, se, sn)

n = length(ve);

% convert lon-lat to xyz
Pxyz = latlon2xyz(dlat, dlon, drad);

% weighting matrix
sig_mag = sqrt(se.^2 + sn.^2);        % magnitude of error
if any(sig_mag == 0)
    Wmat = ones(3,n);
    disp(' NOT using the standard errors in gps2euler.m');
else
    whos s_mag
    Wmat = repmat(1./sig_mag,1,3)';
end

% normalize the input points (METERS)
[Pxyz_norm, Plen] = unit(Pxyz);

% convert to vr-vtheta-vphi
ve  = ve(:);
vn  = vn(:);
vr  = zeros(n,1);
vph = ve;
vth = -vn;

% convert local to global coordinates
Vrtp = [vr vth vph]';
Vxyz = global2local(Vrtp, Pxyz, 0);

% rotated point using small-angle approximation
% (see alternative chunk of code below)
Qxyz = Pxyz + Vxyz;

% normalize the rotated points (METERS)
[Qxyz_norm, Qlen] = unit(Qxyz);

% Pablo Muse: compute the rotation matrix
%whos Qxyz_norm Pxyz_norm Wmat
X = (Wmat .* Qxyz_norm) * (Wmat .* Pxyz_norm)';
[U,S,V] = svd(X);
R = U * diag([1 1 sign(det(U)*det(V)) ]) * V';

% rotation matrix --> euler vector
evec = rotmat2euler(R);   % lat, lon, deg
evec(3) = evec(3) * 1e6;  % deg/yr --> deg/Myr

% convert to wx,wy,wz
evec_xyz = euler_convert(evec,0);
norm(evec_xyz);

% estimated horizontal velocity field
Vrtp = euler2gps(evec_xyz, Pxyz);
ve_est = Vrtp(3,:)' * 1e-3;        % mm/yr --> m/yr
vn_est = -Vrtp(2,:)' * 1e-3;       % mm/yr --> m/yr

%---------------------

if 1==1
    [elat_anti, elon_anti] = antipode(evec(1),evec(2));
    figure; hold on;
    quiver(dlon,dlat,ve_est,vn_est,1,'r');
    plot(evec(2),evec(1),'bo','markersize',10,'markerfacecolor','w');
    plot(elon_anti,elat_anti,'ro','markersize',10,'markerfacecolor','w');
    xlabel('Longitude'); ylabel('Latitude'); title('Rotational field');
    
    figure; hold on;
    quiver(dlon,dlat,ve,vn,1,'b');
    quiver(dlon,dlat,ve_est,vn_est,1,'r');
    xlabel('Longitude'); ylabel('Latitude');
    legend('Input field','Estimated rotational field');
    title(sprintf(' Rotational field with euler vector: (lat = %.4f, lon = %.4f, omega = %.4f deg/Myr)',...
        evec(1),evec(2),evec(3)));

    figure; hold on;
    quiver(dlon,dlat,ve-ve_est,vn-vn_est,'k');
    xlabel('Longitude'); ylabel('Latitude'); title('Residual field');
end

%==============================================================

% % we want to find the DISPLACED point by rotating it from P to Q using the
% % rotation angle specified by the magnitude of the local velocity
% 
% % magnitude of velocity, METERS/yr
% vmag = sqrt( vph.^2 + vth.^2 );
% 
% % rotation angle (RADIANS)
% arc_rot = vmag ./ Plen';
% 
% % compute rotation matrix and apply rotation
% Qxyz = zeros(3,n);
% for ii=1:n
%     % starting point
%     P0 = Pxyz(:,ii);
% 
%     % unit vector for rotation
%     exyz = unit(cross(P0,Vxyz(:,ii)));
%     ex = exyz(1);
%     ey = exyz(2);
%     ez = exyz(3);
% 
%     % rotation angle (radians)
%     omega = arc_rot(ii);
% 
%     % rotation matrix (see euler2rotmat.m)
%     R = zeros(3,3);
%     coso = cos(omega);
%     sino = sin(omega);
%     R(1,1) = ex*ex*(1 - coso) + coso;
%     R(1,2) = ex*ey*(1 - coso) - ez*sino;
%     R(1,3) = ex*ez*(1 - coso) + ey*sino;
%     R(2,1) = ey*ex*(1 - coso) + ez*sino;
%     R(2,2) = ey*ey*(1 - coso) + coso;
%     R(2,3) = ey*ez*(1 - coso) - ex*sino;
%     R(3,1) = ez*ex*(1 - coso) - ey*sino;
%     R(3,2) = ez*ey*(1 - coso) + ex*sino;
%     R(3,3) = ez*ez*(1 - coso) + coso;
% 
%     % apply rotation
%     Qxyz(:,ii) = R * P0;
% end

%==============================================================
