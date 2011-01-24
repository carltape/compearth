%
% function [dlon,dlat,d,dsig,ax0,slabel,ulabel] = get_1D_dataset(ropt,dopt)
% 
% This loads a set of discrete (1D) observations on the sphere
% that we will use to estimate a smooth continuous field.
%
% calls xxx
% called by sphereinterp.m
%

function [dlon,dlat,d,dsig,ax0,slabel,ulabel] = get_1D_dataset(ropt,dopt)

% GEOGRAPHIC REGION 
switch ropt
    case 1, rlabel = 'socal'; ax0 = [-122 -113 30 38]; szone = '11S';
end

%========================================================
% LOAD OBSERVATIONS

dir_data = '/home/carltape/compearth/surfacevel2strain/data/examples/';

switch dopt
    case 1
        dlabel = 'moho03'; 
        ulabel = 'zdep, km';
        [dlon,dlat,d,dsig] = textread([dir_data 'socal_moho_example.dat'],'%f%f%f%f');
end

% subset
inds = getsubset(dlon,dlat,ax0);
dlon = dlon(inds);
dlat = dlat(inds);
d = d(inds);
dsig = dsig(inds);

if isempty(dlon), error('get_1D_dataset.m: zero observations within specified region'); end

% label for files: geographic region, then data set
slabel = [rlabel  '_' dlabel];

ifig = 1;
if ifig==1    
    % plot data and uncertainties
    figure; nr=2; nc=1; msize = 6^2;
    subplot(nr,nc,1); hold on;
    scatter(dlon,dlat,msize,d,'filled');
    %scatter(dlon,dlat,msize,'ko');   % black outline
    colorbar; axis(ax0); grid on;
    xlabel('Longitude'); ylabel('Latitude');
    title(sprintf('%s observations (%i)',slabel,length(d)),'interpreter','none');
    
    subplot(nr,nc,2); hold on;
    scatter(dlon,dlat,msize,dsig,'filled');
    %scatter(dlon,dlat,msize,'ko');   % black outline
    colorbar; axis(ax0); grid on;
    xlabel('Longitude'); ylabel('Latitude');
    title(sprintf('%s uncertainties (%i)',slabel,length(d)),'interpreter','none');

    orient tall, wysiwyg, fontsize(11)
end

%========================================================
