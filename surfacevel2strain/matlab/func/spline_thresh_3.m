%
% function 
% CARL TAPE, 19-Jan-2006
% printed xxx
%
% copied from spline_thresh_2.m on 19-Jan-2006
% See spline_thresh_2.m for options on using convex hull options.
%
% Inputs a matrix and thresholds the columns according to
% whether there are NTRSH entries that exceed the value QTRSH.
% The output is a vector of indices corresponding to the columns
% of A that you want to keep.
%
% calls spline_vals.m
% called by xxx
%

function [ikeep, inum] = spline_thresh_3(spline_tot, qtrsh, ntrsh, dlon, dlat, opts)

% options
ishow = opts{1};

ngrid = length(spline_tot);  % spline_tot is ngrid x 3

k = 0;
ikeep = [];
inum = [];
for ii=1:ngrid   % loop over splines
    
    % spline ii evaluated at all datapoints
    Acol = spline_vals(spline_tot(ii,1), spline_tot(ii,2), spline_tot(ii,3), dlon, dlat, {1});
    
    % datapoint evaluations that are > qtrsh
    itemp = find(Acol > qtrsh);
    nd = length(itemp);
    
    %if ishow==1 disp([ ii nd ]); end
    if nd >= ntrsh
        k = k+1;
        ikeep(k) = ii;
        inum(k)  = length(itemp);
        %if ishow==1, disp([k ii length(itemp) ]); end  
    end
end

ikeep = ikeep(:);
inum = inum(:);

%===================================================