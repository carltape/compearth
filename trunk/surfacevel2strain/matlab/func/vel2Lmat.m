%
% function [D, W, curlv] = vel2Lmat(rtp, vrtp, dvrtpdr, dvrtpdth, dvrtpdph, opts)
% CARL TAPE, 17-May-2006
% printed xxx
%
% Convention is r (up) - theta (south) - phi (east)
%
% This function converts a velocity field in spherical coordinates,
% with gradients, to a velocity gradient tensor, L, where L can be
% decomposed into a symmetric strain tensor, D, and an antisymmetric
% rotation tensor, W.
%    L = D + W, D = D^T, W = -W^T.
%
% Indexing: | 1 2 3 |
%           | 4 5 6 |
%           | 7 8 9 |
%
% See Malvern, p. 671, for spherical formulas.
%
% Ths formulas are singular when th=0 or th=pi.
%
% sopt determines the reduction (if any) of the velocity gradient tensor.
%   sopt = 0, no reduction
%   sopt = 1, reduction based on isotropic linear elasticity
%   sopt = 2, reduction based on viscous rheology
%
% Also computes the curl of the velocity field.
%
% See Latex notes "Supplemental notes fot Tape, Muse, Simons (2008)".
%
% calls xxx
% called by spline_wang_D_figs.m
%

function [D, W, curlv] = vel2Lmat(rtp, vrtp, dvrtpdr, dvrtpdth, dvrtpdph, sopt)

% assign quantities from input parameters
% note the dimensions of the input quantities : n x 3
r       = rtp(:,1);         th      = rtp(:,2);         ph      = rtp(:,3);
vr      = vrtp(:,1);        vth     = vrtp(:,2);        vph     = vrtp(:,3);
dvrdr   = dvrtpdr(:,1);     dvthdr  = dvrtpdr(:,2);     dvphdr  = dvrtpdr(:,3);
dvrdth  = dvrtpdth(:,1);    dvthdth = dvrtpdth(:,2);    dvphdth = dvrtpdth(:,3);
dvrdph  = dvrtpdph(:,1);    dvthdph = dvrtpdph(:,2);    dvphdph = dvrtpdph(:,3);

% cheap way to avoid singularities:
% set all theta values at to at least deps degrees from NP and SP
deps = 1;
ibad = deps*pi/180;
inp = find(th < ibad);    th(inp) = pi/2-ibad;
isp = find(th > pi-ibad); th(isp) = -pi/2+ibad;

num = length(ph);
L = zeros(num,9);
D = zeros(num,9);
W = zeros(num,9);

% Note the singular values for th = 0, ph = 0 .

% L terms WITHOUT radial gradients
% L11
L(:,2) = 1./r .* ( -vth + dvrdth );
L(:,3) = 1./r .* ( -vph + 1./sin(th) .* dvrdph );
% L21
L(:,5) = 1./r .* ( vr + dvthdth );
L(:,6) = 1./r .* ( -vph .* cot(th) + 1./sin(th) .* dvthdph );
% L31
L(:,8) = 1./r .* ( dvphdth );
L(:,9) = 1./r .* ( vr + vth .* cot(th) + 1./sin(th) .* dvphdph );

% L terms WITH radial gradients
if sopt == 0
    disp(' compute L : no reductions');
    L(:,1) = dvrdr;
    L(:,4) = dvthdr;
    L(:,7) = dvphdr;
    
elseif sopt == 1
    disp(' compute L : assume isotropic linear elasticity, Poisson solid, and surface condition');
    F = -1/3;       % Poisson solid (lamda = mu) with F = -lambda / (lambda + 2mu)
    L(:,1) = F * ( L(:,5) + L(:,9) );
    L(:,4) = -L(:,2);
    L(:,7) = -L(:,3);    
    
elseif sopt == 2
    disp(' compute L : assume viscosity and surface condition');
    L(:,1) = zeros(num,1);
    L(:,4) = -L(:,2);
    L(:,7) = -L(:,3);     
    
else
    error('sopt must be 0, 1, or 2');
end

% transpose of L
LT(:,1) = L(:,1);
LT(:,2) = L(:,4);
LT(:,3) = L(:,7);
LT(:,4) = L(:,2);
LT(:,5) = L(:,5);
LT(:,6) = L(:,8);
LT(:,7) = L(:,3);
LT(:,8) = L(:,6);
LT(:,9) = L(:,9);

% symmetric and anti-symmetric parts
D = 0.5*(L + LT);   % six unique elements (1,2,3,5,6,9), D = D^T
W = 0.5*(L - LT);   % three unique elements (2,3,6), W = -W^T
    
% % three components: D_th_th, D_ph_ph, D_th_ph
% % assume zero for all radial direction motion
% D(:,5) = 1./r .* ( dvthdth + vr);
% D(:,9) = 1./r .* ( 1./sin(th) .* dvphdph + vth .* cot(th) );
% D(:,6) = 1./r .* 0.5 .* ( 1./sin(th) .* dvthdph + dvphdth - vph .* cot(th) );
% D(:,8) = D(:,6);

%--------------------------------------------
% CURL of the velocity field (r-theta-phi)

% Malvern p. 670
%curlv = zeros(num,3);
%curlv(:,1) =  1./r .* (vph .* cot(th) + dvphdth - (1./sin(th)) .* dvthdph );
%curlv(:,2) = -1./r .* vph - dvphdr + 1./r .* (1./sin(th)) .* dvrdph ;
%curlv(:,3) =  1./r .* vth + dvthdr - 1./r .* dvrdth;

% Malvern p. 147
curlv = 2 * [-W(:,3) W(:,2) -W(:,1)];

%===================================================
