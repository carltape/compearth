%
% function V_out = global2local(V_in, Pxyz, opts)
% Carl Tape, 03-Sept-2005
%
% This function converts between a global basis and a local basis.
% Pxyz describes points on the sphere (ax,ay,az).
%
% V_in  describes the vector field in (vr,vth,vph) or (vx,vy,vz).
% V_out describes the vector field in (vx,vy,vz) or (vr,vth,vph).
%
% opts:
%   globe2loc==1 then  global -->  local
%   globe2loc==0 then  local  -->  global
%
% calls global2local_rotmat.m
% called by test_global2local.m, plate_model.m
%

function V_out = global2local(V_in, Pxyz, opts)

% ensure that V_in and Pxyz are dimension 3 x n
[m,n1] = size(V_in); if m~=3, V_in = V_in'; end
[m,n2] = size(Pxyz); if m~=3, Pxyz = Pxyz'; end
if n1~=n2, error('Pxyz and V_in must have same dimension'); end

glob2loc = opts(1);

% 'observation' points
[thP, phP] = xyz2tp(Pxyz);

% loop over cartesian components of vectors
V_out = zeros(3,length(thP));
for ii=1:n1
    
    % rotation matrix
    T = global2local_rotmat(thP(ii),phP(ii));
    
    if glob2loc == 1
        % global (x,y,z) to local (r,th,ph)
        V_out(:,ii) = T * V_in(:,ii);
    else
        % local (r,th,ph) to global (x,y,z)
        V_out(:,ii) = inv(T) * V_in(:,ii);
    end
end

%==============================================================
