%
% function dvec = arcdist(slat, slon, rlat, rlon) 
%
% NOTE: JUST USE 'distance' IF YOU HAVE THE MATLAB MAPPING TOOLBOX.
%
% This function determines the arc distance from an
% earthquake source (s) to a seismic receiver (r),
% IGNORING THE OBLATENESS OF THE EARTH (see arcd.m).
%
% INPUT (all in DEGREES):
%    slat = latitude of source
%    slon = longitude of source
%    rlat = latitude of receiver (column vector)
%    rlon = longitude of receiver (column vector)
%
% OUTPUT (in degrees):
%    dvec = arc-distances from (slat, slon) to each (rlat, rlon)
%
% See Matlab function distance.m for geographic coordinates.
%
% calls xxx
% called by many programs
%

function dvec = arcdist(slat, slon, rlat, rlon)

deg = 180/pi;

% make column vectors
rlat = rlat(:); rlon = rlon(:);
slat = slat(:); slon = slon(:);

nrec = length(rlat);
nsrc = length(slat);

if nsrc == 1
    slat = slat * ones(nrec,1);
    slon = slon * ones(nrec,1);
else
    error(' arcdist.m assumes that you have only one source');
end

las = slat / deg;
lar = rlat / deg;
lod = (slon-rlon)/deg ;

dvec = deg * acos( sin(lar).*sin(las) + cos(lar).*cos(las).*cos(lod) );

% Here is the formula when the points are expressed as th, phi
% dvec = acos( cos(th1)*cos(th2) + sin(th1).*sin(th2).*cos(ph1 - ph2) );

%====================================================================
