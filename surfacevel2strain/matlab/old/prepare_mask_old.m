%
% function 
% CARL TAPE, 04-Aug-2007
% printed xxx
%
% Adapted from xy2d.m and c164funA.m on 04-Aug-2007.
%
% 
%
% INPUT:
%   dlon_plot  londitude of plotting point
%   dlat_plot  latitude of plotting point
%   dlon       
%   dlat       
%   D          scalelength for mask
%   Azmin      parameter controlling azimuth threshold
%   Nmin       min number of stations within D degrees and covering Azmin quadrants
%
% OUTPUT:
%
% calls xxx
% called by surfacevel2strain.m
%

function ikeeper = prepare_mask(dlon_plot,dlat_plot,dlon,dlat,D,Nmin,Azmin) 

disp(' computing a mask for plotting...');

nplot = length(dlon_plot);  % numper of plotting point
ndata = length(dlon);       % number of datapoints

% if the desired number of points is less than the desired number of
% quadrants covered, then this is impossible
if Nmin < Azmin, Nmin = Azmin; end

% number of azimuth bins
az_bins = 4;

bikeeper = zeros(nplot,1);
for ii = 1:nplot
    if mod(ii,100) == 0, disp(sprintf('%i out of %i',ii,nplot)); end
    
    % target plotting point
    lon0 = dlon_plot(ii);
    lat0 = dlat_plot(ii);
    %figure; plot(dlon,dlat,'b.',lon0,lat0,'rx'); axis equal

    % compute angular distance and azimuth to each datapoint
    [dist0,az0] = distance(lat0,lon0,dlat,dlon);

    % number of datapoints within the critical angular distance
    igood = find(dist0 <= D);
    ngood = length(igood);
    
    % bin these datapoints into four azimuthal quadrants
    Nbin = zeros(5,1);
    if length(igood) > Nmin
        Nbin = histc(az0(igood),linspace(0,360,az_bins+1));
        if sum(sign(Nbin)) >= Azmin
            bikeeper(ii) = 1;
        end
    end
end

ikeeper = find( bikeeper == 1);

%================================================================
