function utm2ll_all_ellipsoids(xi,yi,szone,itype)
%UTM2LL_ALL_ELLIPSOIDS check all available ellipsoids for a set of points
% This is a function that may be useful in comparing the calculations from
% Matlab's mapping toolbox with calculations from another code (or source).
%
% EXAMPLES: utm2ll_all_ellipsoids(178.40,-38.92,'60H',0)
%           utm2ll_all_ellipsoids(621367.42,5691169.40,'60H',1)
%

% see the reference page for the alamanac function
etags = {'everest','bessel','airy','clarke66','clarke80','international',...
'krasovsky','wgs60','iau65','wgs66','iau68','wgs72','grs80','wgs84'};
ne = length(etags);

if itype==0, stfmt = '%.2f %.2f'; else stfmt = '%.10f %.10f'; end

n = length(xi);

for ii=1:n
   fprintf('Point %i/%i at (%.6f, %.6f)\n',ii,n,xi(ii),yi(ii));
   for kk=1:ne
        ellipsoid = almanac('earth',etags{kk},'meters');
        [ex,ey] = utm2ll(xi(ii),yi(ii),szone,itype,ellipsoid);
        fprintf(['%14s ellipsoid (utm zone %s) --> ' stfmt '\n'],etags{kk},szone,ex,ey);
   end
end
