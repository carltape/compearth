%
% function [G, Gdph, Gdth] = dogsph_vals_mat(spline_tot, dlon, dlat)
% Pablo Muse, 04-April-2007
%
% Given a set of spherical splines and datapoints, this function returns
% the design matrix for estimating functions on a sphere.
%
% spline_tot : N x 3 array with columns lon-center, lat-center, q-grid
% dlon, dlat : locations of datapoints for the inverse problem
% ndim       : number of components for each observation
%              ndim = 1 for scalar field
%              ndim = 3 for 3-component vector field
%
% NOTE: This function could also be written to take in a model vector for a
% scalar field, and then output the VECTORS g, dgdhh, and dgdth.
%
% calls dogsph_vals.m
% called by surfacevel2strain.m
%

%function [G, Gdph, Gdth] = spline_vals_mat(spline_tot, dlon, dlat, ndim)
function [G, Gdph, Gdth, Gdphdph, Gdthdth, Gdthdph] = dogsph_vals_mat(spline_tot, dlon, dlat, opts)

ndata = length(dlon);
ngrid = length(spline_tot);

G = zeros(ndata, ngrid);
Gdph = zeros(ndata, ngrid);
Gdth = zeros(ndata, ngrid);

if nargout == 6
    Gdphdph = zeros(ndata, ngrid);
    Gdthdth = zeros(ndata, ngrid);
    Gdthdph = zeros(ndata, ngrid);
end
% fill each column of G with a basis function evaluated at all the datapoints
for jj=1:ngrid
    %disp(jj);
    ff = dogsph_vals(spline_tot(jj,1), spline_tot(jj,2), spline_tot(jj,3), dlon, dlat, opts);
    G(:,jj) = ff(:,1);
    Gdph(:,jj) = ff(:,2);
    Gdth(:,jj) = ff(:,3);
    if nargout == 6
        Gdphdph(:,jj) = ff(:,4);
        Gdthdth(:,jj) = ff(:,5); 
        Gdthdph(:,jj) = ff(:,6); 
    end   
end

% if ndim == 1
%     G = G0;
%     Gdph = G0dph;
%     Gdth = G0dth;
% else
%     % fill cell arrays
%     for jj=1:ngrid
%         for ii=1:ndata
%             Gcell{ii,jj} = G0(ii,jj) * eye(ndim);
%             Gcelldph{ii,jj} = G0dph(ii,jj) * eye(ndim);
%             Gcelldth{ii,jj} = G0dth(ii,jj) * eye(ndim);
%         end
%     end
%     
%     % convert cell arrays to numeric arrays
%     G = zeros(ndata*ndim, ngrid*ndim);
%     Gdph = zeros(ndata*ndim, ngrid*ndim);
%     Gdth = zeros(ndata*ndim, ngrid*ndim);
%     for ii = 1:ndata
%         Gii = [];
%         for jj = 1:ngrid
%             Gii = [Gii Gcell{ii,jj}];
%         end
%         irows = ndim*ii-(ndim-1) : ndim*ii;
%         G(irows,:) = Gii;
%     end
% end

%=====================================================================
