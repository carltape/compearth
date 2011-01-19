%
% function Elatlon = poles2euler(Plat1,Plon1,Plat2,Plon2)
% Carl Tape, 20-June-2008
% printed xxxx
%
% This function takes two latitude-longitude points in degress,
% and computes the pole of rotation.
%
% Example (see below).
%
% calls xxx
% called by xxx
%

function Elatlon = poles2euler(Plat1,Plon1,Plat2,Plon2)

[Exyz,Elat,Elon] = latlons2pole(Plat1,Plon1,Plat2,Plon2);

[Pdist, az] = distance(Plat1,Plon1,Plat2,Plon2);

%    evec(1) = latitude (deg) of euler pole
%    evec(2) = longitude (deg) of euler pole
%    evec(3) = rotation angle (deg)
Elatlon = [Elat ; Elon ; Pdist];

%--------------

if 0==1
    % specify the starting pole P1
    Plat1 = 10; Plon1 = 80;
    
    % specify the finishing pole P2 by getting the pole of a great circle
    lat1 = 36; lon1 = -120; lat2 = 34; lon2 = -117;
    [Pxyz,Plat2,Plon2] = latlons2pole(lat1,lon1,lat2,lon2);
    
    % get the pole describing P1 to P2
    Elatlon = poles2euler(Plat1,Plon1,Plat2,Plon2);
    
    % CHECK: apply the rotation to P1, then compare P2 and the rotated P1
    [lat_rot, lon_rot, R] = euler_rot_tec(Plat1,Plon1,Elatlon);
    Plat2, Plon2
    lat_rot, lon_rot
end

%===========================================================
