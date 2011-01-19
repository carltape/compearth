%
% function [D, W, curlv] = vel2rot(rtp, vrtp, dvrtpdr, dvrtpdth, dvrtpdph, opts)
% CARL TAPE, 17-May-2006
% printed xxx
%
% Copied from vel2Lmat.m on 21-June-2009.
%
% calls xxx
% called by spline_wang_D_figs.m
%

function [ws_r,ws_th,ws_ph,wt_r,wt_th,wt_ph] = vel2rot(rtp, vrtp, dvrtpdr, dvrtpdth, dvrtpdph)

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
ws_r  = zeros(num,1);
ws_th = zeros(num,1);
ws_ph = zeros(num,1);
wt_r  = zeros(num,1);
wt_th = zeros(num,1);
wt_ph = zeros(num,1);

% rotation components associated with vtheta and vphi
ws_r  = 1./(2*r) .* ( vph.*cot(th) - 1./sin(th).*dvthdph + dvphdth );
ws_th = 1./r     .* ( -vph );
ws_ph = 1./r     .* ( vth );

% rotation components associated with vr
wt_r  = zeros(num,1);
wt_th = 1./r     .* ( 1./sin(th) .* dvrdph );
wt_ph = 1./r     .* ( -dvrdth );

%===================================================
