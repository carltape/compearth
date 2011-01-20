%
% function strike_slip_sphere.m
% Carl Tape, 24-July-2008
%
% This generates the displacement field for an equatorial great-circle
% strike-slip fault by approximating the Okada solution.
%
% INPUT
%   dlat,dlon       observation pointsm degrees
%   dep             locking depth, meters
%   amp             amplitude of slip rate, meters/yr
%   lat1,lon1       starting point of fault
%   lat2,lon2       finishing point of fault
%
% OUTPUT
%   Vtheta,Vphi     local velocity, meters/yr
%   dmin            distance to pole of great circle, degrees
%   Glat,Glon       lat,lon of great circle
%   R
%   evecA,evecB     euler vectors of two shells (wx,wy,wx)
%
% EXAMPLES (see below).
%
% calls latlons2pole.m, euler_rot_tec.m, poles2euler.m
% called by xxx
%

function [Vtheta_mod,Vphi_mod,dmin,Glat,Glon,evecA,evecB] = ...
    strike_slip_sphere(dlat,dlon,dep,amp,lat1,lon1,lat2,lon2)

rad = 6371*1e3;
deg = 180/pi;
nobs = length(dlon);
ax0 = [min(dlon) max(dlon) min(dlat) max(dlat)];

disp('computing velocity field for a locked great-circle strike-slip fault...');

% get input points describing the strike-slip fault
n = 1000;
gc_lon = linspace(-180,180,n)';     % note the range in COLATITUDE
gc_lat = zeros(n,1);

% specify the starting pole P1 of the great-circle fault (North Pole)
Plat1 = 90; Plon1 = 0;
    
% compute the finishing pole P2 of the fault by getting the pole of the great circle
[Pxyz2,Plat2,Plon2] = latlons2pole(lat1,lon1,lat2,lon2);

% get the Euler vector describing P1 to P2
Elatlon = poles2euler(Plat1,Plon1,Plat2,Plon2);

% apply rotation to the great circle -- Glat and Glon now describe
% the great-circle (or plane) of the desired strike-slip fault
[Glat, Glon, R] = euler_rot_tec(gc_lat,gc_lon,Elatlon);

%figure; plot(gc_lon,gc_lat,'b.',Glon,Glat,'r.');
%legend('equator points','rotated great circle');
%xlabel('Longitude, deg'); ylabel('Latitude, deg');

%----------------------
% PART 1: consider slip rate

% compute distance from observation point to the pole of the great-circle fault
dpole = zeros(nobs,1);
dpole = distance(repmat(Plat2,nobs,1), repmat(Plon2,nobs,1), dlat, dlon);
dmin = 90 - dpole;

% for points <90 deg from the pole, give them positive rotation
% for points >90 deg from the pole, give them negative rotation
iA = find( dpole <= 90 );
iB = find( dpole > 90 );

% get xyz of observation points
Pxyz = latlon2xyz(dlat,dlon,rad);

% angular distance to move in one MILLION years, in deg/Myr
% NOTE: input slip rate is split evenly by the FACTOR OF TWO
amp_theta = amp / rad * deg * 1e6 / 2;

% compute euler vectors
evecA0 = [Plat2 Plon2 amp_theta];
evecB0 = [Plat2 Plon2 -amp_theta];
evecA = euler_convert(evecA0,2);
evecB = euler_convert(evecB0,2);

% compute local velocities: mm/yr --> m/yr
VrtpA = 1e-3 * euler2gps(evecA, Pxyz(:,iA));
VrtpB = 1e-3 * euler2gps(evecB, Pxyz(:,iB));

Vtheta = zeros(nobs,1);
Vtheta(iA) = VrtpA(2,:)';
Vtheta(iB) = VrtpB(2,:)';

Vphi = zeros(nobs,1);
Vphi(iA) = VrtpA(3,:)';
Vphi(iB) = VrtpB(3,:)';

figure; hold on;
quiver(dlon,dlat,Vphi,-Vtheta,1,'b');

%----------------------
% PART 2: consider locking depth -- this will modified the magnitude of the
% vector field

if dep == 0
    Vmod = ones(nobs,1);                % no locking (creeping fault!)
else
    % dpole is [0,180]
    % Vmod is [-1,1]
    % abs value is because the input field already has the negative component
    Vmod = (2/pi) * abs( atan2(rad*(dpole-90)/deg, dep) );
end

Vmag = sqrt( Vtheta.^2 + Vphi.^2 );
Vunit = [Vtheta  Vphi] ./ repmat(Vmag,1,2);
Vmag_mod = Vunit .* repmat(Vmag .* Vmod,1,2);
Vtheta_mod = Vmag_mod(:,1);
Vphi_mod = Vmag_mod(:,2);

quiver(dlon,dlat,Vphi_mod,-Vtheta_mod,1,'r');
legend('unlocked',['locked at ' num2str(dep/1e3) ' km depth']);
plot(Glon,Glat,'k.','markersize',4); axis equal; axis(ax0);
xlabel('Latitude'); ylabel('Longitude');

figure; quiver(dlon,dlat,Vphi-Vphi_mod,-Vtheta+Vtheta_mod,'k');
title('Residual between locked and unlocked velocity fields');
xlabel('Latitude'); ylabel('Longitude'); axis(ax0);

%----------------------

% % compute shortest distance from each observation point to the great-circle fault
% % NOTE: probably this should be done ANALYTICALLY, but it will bo okay if
% %       the parameterization of the input fault is dense enough (see n above)
% dmin = zeros(nobs,1);
% azi = zeros(nobs,1);
% for ii = 1:nobs
%     [dall,azall] = distance(repmat(dlat(ii),n,1), repmat(dlon(ii),n,1), Glat, Glon);
%     [junk,imin] = min(dall);
%     azi(ii) = azall(imin);
%     dmin(ii) = dall(imin);
%     %disp([ii dlat(ii) dlon(ii) azi(ii) dmin(ii)]);
% end
% 
% % compute delta angles, which are in COLATITUDE
% iSW = (azi > 180);
% dmin(iSW) = -dmin(iSW);
% deltas = 90 - dmin;
% 
% % compute the displacement field as if the great circle were the equator
% [Utheta0,Uphi0] = strike_slip_equator(rad,dep,amp,deltas);
% 
% % convert displacement field to global coordinates
% Vrtp = [zeros(nobs,1) Utheta0 Uphi0]';
% Pxyz = latlon2xyz(dlat,dlon,rad);
% Vxyz = global2local(Vrtp,Pxyz,0);
% 
% % rotate global velocity field from equator reference to fault reference
% Vrot_xyz = R*Vxyz;
% 
% % convert to local coordinates
% Vrot_rtp = global2local(Vrot_xyz,Pxyz,1);
% Vtheta = Vrot_rtp(2,:)';
% Vphi = Vrot_rtp(3,:)';
% 
% figure; hold on;
% quiver(dlon,dlat,Vphi,-Vtheta); axis equal
% plot(Glon,Glat,'r.');

%----------------------

if 0==1
    %----------------------------
    % EXAMPLE 1: socal points
    
    % specify input points
    clear, close all, clc
    ax0 = [-122 -113 30 38];
    lonmin = ax0(1); lonmax = ax0(2);
    latmin = ax0(3); latmax = ax0(4);
    [dlon,dlat] = gridvec(lonmin,lonmax,30,latmin,latmax);
    %figure; plot(dlon,dlat,'.');
    
    dep = 15*1e3;       % locking depth, m
    slip = 0.04;        % FULL slip rate, m/yr
    lat1 = 36; lon1 = -120; lat2 = 34; lon2 = -117;
    [Vtheta,Vphi,dmin,evecA,evecB,Glat,Glon] = ...
        strike_slip_sphere(dlat,dlon,dep,slip,lat1,lon1,lat2,lon2);
    
    % distance to the fault, in meters
    xdat = 6371*1e3*dmin*pi/180;
    
    Vmag = sqrt( Vtheta.^2 + Vphi.^2 );
    figure; plot(xdat, Vmag*1e3,'b.'); xlabel(' observation index'); ylabel('speed, mm/yr');
    hold on;
    
    % check the magnitudes computed from the old 2D code
    [dh,dv] = fault_2D(xdat,slip,Inf,90,dep,1);
    plot(xdat,dh*1e3,'r.');
    
    %----------------------------
    % EXAMPLE 2: global points
    
    % specify input points
    clear, close all, clc
    ax0 = [-180 180 -90 90];
    lonmin = ax0(1); lonmax = ax0(2);
    latmin = ax0(3); latmax = ax0(4);
    [dlon,dlat] = gridvec(lonmin,lonmax,100,latmin,latmax);
    
    dep = 300*1e3;       % locking depth, m
    amp = 0.04;        % FULL slip rate, m/yr
    lat1 = 36; lon1 = -120; lat2 = 34; lon2 = -117;
    [Vtheta,Vphi,dmin,evecA,evecB,Glat,Glon] = ...
        strike_slip_sphere(dlat,dlon,dep,amp,lat1,lon1,lat2,lon2);
end

%==============================================================
