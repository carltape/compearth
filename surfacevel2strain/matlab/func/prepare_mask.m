%
% function prepare_mask.m
% Carl Tape, 13-Nov-2008
%
% This function uses a very specific, somewhat automated procedure to
% compute a mask for obscuring estimated points with large posterior
% uncertainties.
%
% calls xxx
% called by surfacevel2strain.m
%

function [ikeeper,sigval,sigval_cum] = prepare_mask(dlon_plot,dlat_plot,G_plot,Cm_s,Cm_e,Pcum) 

disp(' computing a mask for plotting...');

nplot = length(dlon_plot);  % numper of plotting point

% compute diagonal of posterior data covariance matrix
Cd_s_diag_plot = diag(G_plot*Cm_s*G_plot');   % nplot x 1
Cd_e_diag_plot = diag(G_plot*Cm_e*G_plot');   % nplot x 1
sn_post_plot = sqrt(Cd_s_diag_plot);          % nplot x 1
se_post_plot = sqrt(Cd_e_diag_plot);          % nplot x 1

% compute a scalar function to represent the estimated error
sigval = log( sn_post_plot .* se_post_plot * 1e6 );  % ln of the area of ellipse

% use sigval to select the good points, i.e., the points NOT to mask out
[~, isigsort] = sort(sigval);
sigval_shift = sigval - min(sigval);
sigval_norm = sigval_shift / sum(sigval_shift);
sigval_cum = cumsum(sort(sigval_norm) );
ikeeper = isigsort(find(sigval_cum < Pcum));

%------------------------------------------------------------
% old version

% function ikeeper = prepare_mask(dlon_plot,dlat_plot,dlon,dlat,D,Nmin,Azmin) 
% 
% disp(' computing a mask for plotting...');
% 
% nplot = length(dlon_plot);  % numper of plotting point
% ndata = length(dlon);       % number of datapoints
% 
% % if the desired number of points is less than the desired number of
% % quadrants covered, then this is impossible
% if Nmin < Azmin, Nmin = Azmin; end
% 
% % number of azimuth bins
% az_bins = 4;
% 
% bikeeper = zeros(nplot,1);
% for ii = 1:nplot
%     if mod(ii,100) == 0, disp(sprintf('%i out of %i',ii,nplot)); end
%     
%     % target plotting point
%     lon0 = dlon_plot(ii);
%     lat0 = dlat_plot(ii);
%     %figure; plot(dlon,dlat,'b.',lon0,lat0,'rx'); axis equal
% 
%     % compute angular distance and azimuth to each datapoint
%     [dist0,az0] = distance(lat0,lon0,dlat,dlon);
% 
%     % number of datapoints within the critical angular distance
%     igood = find(dist0 <= D);
%     ngood = length(igood);
%     
%     % bin these datapoints into four azimuthal quadrants
%     Nbin = zeros(5,1);
%     if length(igood) > Nmin
%         Nbin = histc(az0(igood),linspace(0,360,az_bins+1));
%         if sum(sign(Nbin)) >= Azmin
%             bikeeper(ii) = 1;
%         end
%     end
% end
% 
% ikeeper = find( bikeeper == 1);

%================================================================
