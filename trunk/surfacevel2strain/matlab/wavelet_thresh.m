% function [ikeep, inum] = wavelet_thresh(gridpoints,
% qtrsh_cos_dist, ntrsh, dlon, dlat)
% dlon, dlat)
%
% CARL TAPE's code SPLINE_THRESH_3.M modified by Pablo Muse, 
% to be used with spherical wavelets instead of spherical splines. 
%
% PM, 28-Aug-2007
%
%
% Inputs a matrix and thresholds the columns according to
% whether there are NTRSH entries that exceed the value QTRSH.
% The output is a vector of indices corresponding to the columns
% of A that you want to keep.
%
% calls dogsph_vals.m
% called by xxx
%

function [ikeep, inum] = wavelet_thresh(gridpoints, qtrsh_cos_dist, ntrsh, dlon, dlat)

% options
%ishow = opts{1};

ngrid = length(gridpoints);  % spline_tot is ngrid x 3

k = 0;
ikeep = [];
inum = [];
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

is_in_support = zeros(size(dlon));

for ii=1:ngrid   % loop gridpoints
    %if (gridpoints(ii,3) >= 1)    % commented out on 19-June-2008
        is_in_support = (dataptsxyz*(gridptsxyz(ii,:))' > qtrsh_cos_dist(gridpoints(ii,3)+1));
        nb_in = sum(is_in_support);
        if (nb_in >= ntrsh)
            k = k+1;
            ikeep(k) = ii;
            inum(k)  = nb_in;
        end
    %end
end

ikeep = ikeep(:);
inum = inum(:);

%===================================================