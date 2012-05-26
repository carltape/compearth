function Mout = global2local_mat(Min,lat,lon,iglob2loc)
%GLOBAL2LOCAL_MAT converts a set of matrices between global and local bases
%
%
% INPUT
%   Min         3 x 3 x n set of matrices
%   lat,lon     input points on the sphere
%   iglob2loc   =1 for (x,y,z) --> (r,th,ph)
%               =0 for (r,th,ph) --> (x,y,z)
% OUTPUT
%   Mout        3 x 3 x n set of matrices in new basis
%
% The vector version of this function is global2local.m.
%
% calls global2local_rotmat.m
% called by global2local_MT.m
%

% ensure that Min and Pxyz have the same number of entries
[a,b,n1] = size(Min);
if any([a b]~=3), error('M should be 3 x 3 x n'); end
if n1 ~= length(lat), error('M should have same length as lat, lon'); end

% compute theta, phi for each point
deg = 180/pi;
thP = (90 - lat)/deg;
phP = lon/deg;

% loop over input points
Mout = 0*Min;
for ii=1:n1
    
    % rotation matrix
    % note 1: this will DIFFER for each input point
    % note 2: this  will be orthogonal: T T^t = I and T^t = T^-1
    T = global2local_rotmat(thP(ii),phP(ii));
    
    % input matrix
    M0 = Min(:,:,ii);
    
    if iglob2loc == 1
        % global (x,y,z) to local (r,th,ph)
        %Mout(:,:,ii) = T * M0 * inv(T);
        Mout(:,:,ii) = T * M0 * T';
    else
        % local (r,th,ph) to global (x,y,z)
        %Mout(:,:,ii) = inv(T) * M0 * T;
        Mout(:,:,ii) = T' * M0 * T;
    end
end

%==============================================================
