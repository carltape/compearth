function dvec = sdistance(lat1,lon1,lat2,lon2)
% SDISTANCE arc distance between sets of points on the sphere
%
% OUTPUT (in degrees):
%    dvec = arc-distances from (lat1, lon1) to each (lat2, lon2)
%
% NOTE: JUST USE 'distance' IF YOU HAVE THE MATLAB MAPPING TOOLBOX;
% THIS WILL ALSO PROVIDE THE AZIMUTH AND PROVIDE OPTIONS FOR GEOGRAPHIC
% CORRDINATES.
%

deg = 180/pi;

% make column vectors
lat1 = lat1(:); lon1 = lon1(:);
lat2 = lat2(:); lon2 = lon2(:);

if length(unique([length(lat1) length(lon1) length(lat2) length(lon2)]))~=1
    whos lat1 lon1 lat2 lon2
    error('input vectors must all have the same length');
end

las = lat2 / deg;
lar = lat1 / deg;
lod = (lon2-lon1)/deg ;
dvec = deg * acos( sin(lar).*sin(las) + cos(lar).*cos(las).*cos(lod) );

% Here is the formula when the points are expressed as th, phi
% dvec = acos( cos(th1)*cos(th2) + sin(th1).*sin(th2).*cos(ph1 - ph2) );

%==========================================================================

if 0==1
   n = 10;
   lat1 = -90 + 180*rand(n,1); 
   lon1 = -180 + 360*rand(n,1); 
   lat2 = -90 + 180*rand(n,1); 
   lon2 = -180 + 360*rand(n,1); 
   dvec = sdistance(lat1,lon1,lat2,lon2);
   dmap = distance(lat1,lon1,lat2,lon2);    % matlab mapping toolbox
   [dvec dmap]
end

%==========================================================================