% function [ikeep, inum] = wavelet_min_scale(gridpoints, qtrsh_cos_dist, ntrsh, dlon, dlat)
%
% INPUT
%    gridpoints      ngrid x 3 matrix of (glon, glat, qscale)
%    qtrsh_cos_dist  cosine of angular footprint of each q-scale spherical wavelet
%    dlon, dlat      ndata x 1 vectors of input points
%
% OUTPUT
%    qmin_scale      ndata x 1 vector of minimum q footprint touching each point
%
% Modified from wavelet_thresh.m on 14-June-2009
%
% calls xxx
% called by xxx
%

function qmin_scale = wavelet_min_scale(gridpoints, qtrsh_cos_dist, dlon, dlat)

ndata = length(dlon);
ngrid = length(gridpoints);  % spline_tot is ngrid x 3
rad = pi/180;

% data points
th = (90 - dlat).*rad;
ph = dlon.*rad;
sin_th = sin(th);
cos_th = cos(th);
sin_ph = sin(ph);
cos_ph = cos(ph);

% grid points
th_gridpts = (90 - gridpoints(:,2))*rad;
ph_gridpts = gridpoints(:,1)*rad;
gridpts_sin_th = sin(th_gridpts);
gridpts_cos_th = cos(th_gridpts);
gridpts_sin_ph = sin(ph_gridpts);
gridpts_cos_ph = cos(ph_gridpts);

% convert to xyz
dataptsxyz = [sin_th.*cos_ph sin_th.*sin_ph cos_th];
gridptsxyz = [gridpts_sin_th.*gridpts_cos_ph gridpts_sin_th.*gridpts_sin_ph gridpts_cos_th];

% support of each spherical wavelet
qsupp_cos = qtrsh_cos_dist(gridpoints(:,3)+1);

qmin_scale = NaN*ones(ndata,1);
for ii=1:ndata   % loop datapoints
    is_in_support = zeros(ngrid,1);
    is_in_support = (gridptsxyz*(dataptsxyz(ii,:))' > qsupp_cos);
    qmax = max( gridpoints(is_in_support,3) );
    qmin_scale(ii) = qmax;
end

%===================================================