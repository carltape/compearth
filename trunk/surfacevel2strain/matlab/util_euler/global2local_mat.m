%
% function Mout = global2local_mat(Min,Pxyz,opts)
% Carl Tape, 15-March-2011
%
% This function converts a set of SYMMETRIC MATRICES between
% a global basis and a local basis.
% Pxyz describes points on the sphere (ax,ay,az).
%
% Min  describes the matrix in (r,th,ph) or (x,y,z).
% Mout describes the matrix in (x,y,z) or (r,th,ph).
%
% opts:
%   globe2loc==1 then  global -->  local
%   globe2loc==0 then  local  -->  global
%
% FUTURE WORK: Why not have it work for non-symmetric input matrices (9 x 9)?
%
% The vector version of this function is global2local.m.
%
% calls global2local_rotmat.m
% called by rotate_points_MT.m
%

function [Mout,Tall] = global2local_mat(Min,Pxyz,opts)

% check that M is 6 x n
[Min,n1] = Mdim(Min);

% transform to 3 x 3 x n
Min = Mvec2Mmat(Min,1);

% ensure that Min and Pxyz are dimension 3 x n
[m,n2] = size(Pxyz); if m~=3, Pxyz = Pxyz'; end
if n1~=n2, error('Pxyz and Min must have same lengths'); end

glob2loc = opts(1);

% 'observation' points
[thP, phP] = xyz2tp(Pxyz);

% loop over cartesian components of vectors
Mout = 0*Min;
Tall = 0*Min;
for ii=1:n1
    
    % rotation matrix
    % note 1: this will DIFFER for each input point
    % note 2: this  will be orthogonal, T = T^t = T^-1
    T = global2local_rotmat(thP(ii),phP(ii));
    Tall(:,:,ii) = T;
    
    % input matrix
    M0 = Min(:,:,ii);
    
    if glob2loc == 1
        % global (x,y,z) to local (r,th,ph)
        %Mout(:,:,ii) = T * M0 * inv(T);
        Mout(:,:,ii) = T * M0 * T';
    else
        % local (r,th,ph) to global (x,y,z)
        %Mout(:,:,ii) = inv(T) * M0 * T;
        Mout(:,:,ii) = T' * M0 * T;
    end
end

% transform to 6 x n
Mout = Mvec2Mmat(Mout,0);

%==============================================================
