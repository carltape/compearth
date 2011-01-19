%
% function 
% CARL TAPE, 13-Nov-2008
% printed xxx
%
%
% 
%
% INPUT:
%   dlon_plot  londitude of plotting point
%   dlat_plot  latitude of plotting point
%
% OUTPUT:
%
% calls xxx
% called by surfacevel2strain.m
%

function [ikeeper,sigval,sigval_cum] = prepare_mask(dlon_plot,dlat_plot,G_plot,Cm_s,Cm_e,Pcum) 

disp(' computing a mask for plotting...');

nplot = length(dlon_plot);  % numper of plotting point

% compute diagonal of posterior data covariance matrix
Cd_s_diag_plot = diag(G_plot*Cm_s*G_plot');   % nplot x 1
Cd_e_diag_plot = diag(G_plot*Cm_e*G_plot');   % nplot x 1
sn_post_plot = sqrt(Cd_s_diag_plot);          % nplot x 1
se_post_plot = sqrt(Cd_e_diag_plot);          % nplot x 1

% compute a scalar function to represent the estimated error
sigval = log( sn_post_plot .* se_post_plot * 1e6 );  % ln of the area of ellipse

% use sigval to select the good points, i.e., the points NOT to mask out
[junk, isigsort] = sort(sigval);
sigval_shift = sigval - min(sigval);
sigval_norm = sigval_shift / sum(sigval_shift);
sigval_cum = cumsum(sort(sigval_norm) );
ikeeper = isigsort(find(sigval_cum < Pcum));



%================================================================
