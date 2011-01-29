%
% function [spline_tot, G, Gdph, Gdth] = spline_vals_mat(spline_tot, dlon, dlat)
% Carl Tape, 19-Sept-2006
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
% calls spline_vals.m
% called by test_platemodel2strain.m
%

%function [G, Gdph, Gdth] = spline_vals_mat(spline_tot, dlon, dlat, ndim)
function [G, Gdph, Gdth] = spline_vals_mat(spline_tot, dlon, dlat)

ndata = length(dlon);
ngrid = length(spline_tot);

   G0 = zeros(ndata, ngrid);
G0dph = zeros(ndata, ngrid);
G0dth = zeros(ndata, ngrid);

% fill each column of G with a basis function evaluated at all the datapoints
for jj=1:ngrid
    %disp(jj);
    ff = spline_vals(spline_tot(jj,1), spline_tot(jj,2), spline_tot(jj,3), dlon, dlat, {3});
    
       G0(:,jj) = ff(:,1);
    G0dph(:,jj) = ff(:,2);
    G0dth(:,jj) = ff(:,3);
end

G = G0;
Gdph = G0dph;
Gdth = G0dth;

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
