%
% forward_mogi.m
% Carl Tape, 08-Aug-2007
%
% Mogi source used in Tape et al. (GJI 2009).
%
% calls utm2ll.m
% called by socal_gps_syn.m
%

function [Ux, Uy, Uz] = forward_mogi(m,lon,lat,szone)

% model vector describing center of the volcanic source
lon0   = m(1);
lat0   = m(2);
depth0 = m(3);
c      = m(4);

% convert lon-lat to utm (all points are in meters)
[x0, y0] = utm2ll(lon0,lat0,szone,0);
[x, y]   = utm2ll(lon,lat,szone,0);
    
dx = x-x0;
dy = y-y0;
dz = depth0;

r  = sqrt(dx.^2 + dy.^2 + dz.^2);

Ux = c * dx ./ (r.^3);
Uy = c * dy ./ (r.^3);
Uz = c * dz ./ (r.^3);

Ux = Ux(:);
Uy = Uy(:);
Uz = Uz(:);

%================================================================
