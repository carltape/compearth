%
% socal_gps_split.m
% Carl Tape, 16-Jan-2008
%
%
% calls xxx
% called by xxx
%

function [iwest, ieast] = socal_gps_split(ax0,lon,lat)

% constants
deg = 180/pi;
earthr = 6371*1e3;
ifig = 1;

lonmin = ax0(1);
lonmax = ax0(2);
latmin = ax0(3);
latmax = ax0(4);

% USER: load the SAF for plotting: latsaf, lonsaf, xsaf, ysaf
% NOTE: San Andreas fault is indexed from NORTH to SOUTH
gdir = '/home/carltape/compearth/surfacevel2strain/gmt/input/';
%gdir = '../gmt/input/';
if exist(gdir)==7
    [lonsaf,latsaf,xsay,ysaf] = textread([gdir 'safdata2.dat'],'%f%f%f%f');
    nsaf = length(lonsaf);
else
    error(['gdir does not exist: ' gdir]);
end
    
% segment of the SAF that lies within the box
%isub0 = find(and( latsaf >= latmin, latsaf <= latmax ));
[isaf_in, isaf_on] = inpolygon(lonsaf,latsaf,ax0([1 2 2 1 1]),ax0([3 3 4 4 3]));
isub0 = find(isaf_in == 1);
%isub0 = getsubset(lonsaf,latsaf,ax0);
isub = [isub0(1)-1 ; isub0 ; isub0(end)+1];   % add a point on both ends
blon = lonsaf(isub);
blat = latsaf(isub);

% Pacific plate "boundary", starting from SW corner going clockwise
pa_lon = [lonmin ; blon ; lonmin];
pa_lat = [latmin ; blat ; latmin];

% North America plate "boundary", starting from NE corner going counter-clockwise
na_lon = [lonmax ; lonmin ; blon ; lonmax ; lonmax];
na_lat = [latmax ; latmax ; blat ; latmin ; latmax];

%-------------------------------------------------------

[ina_in, ina_on] = inpolygon(lon,lat,na_lon,na_lat);
[ipa_in, ipa_on] = inpolygon(lon,lat,pa_lon,pa_lat);
[igps_in, igps_on] = inpolygon(lon,lat,ax0([1 2 2 1 1]),ax0([3 3 4 4 3]));
ieast = find(ina_in == 1);
iwest = find(ipa_in == 1);
igps  = find(igps_in == 1);

if length(iwest)+length(ieast) ~= length(igps), error('mismatch of indexing'); end

if ifig == 1
    figure; hold on;
    plot(lonsaf, latsaf, 'bo')
    plot(lonsaf(isub), latsaf(isub), 'ro')
    plot(ax0([1 2 2 1 1]),ax0([4 4 3 3 4]),'k')
    plot(na_lon, na_lat, 'b.--')
    plot(pa_lon, pa_lat, 'r.--')

    figure; hold on; plot(na_lon, na_lat, 'b.--'); plot(lon(ina_in), lat(ina_in), 'b+');
    figure; hold on; plot(pa_lon, pa_lat, 'r.--'); plot(lon(ipa_in), lat(ipa_in), 'r+');
end

%========================================================
