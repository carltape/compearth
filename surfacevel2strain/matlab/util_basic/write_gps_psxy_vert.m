%
% function write_gps_psxy_vert.m
%
% This file writes the vertical components of a GPS velocity field to a
% file for plotting in GMT using psxy.
%
% calls xxx
% called by xxx
%

function write_gps_psxy_vert(filetag,lon,lat,vu,su)

filename = [filetag '_psxy_vert.dat'];
nstation = length(lon);

% sort in by decreasing dot size, so that small dots are not buried by
% large ones in the psxy plot
[~,isort] = sort(abs(vu),'descend');

% COLUMNS : lon lat vu su
fid = fopen(filename,'w');
for ii = 1:nstation
    jj = isort(ii);
%    fprintf(fid,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n', lon(ii),lat(ii),vu(ii),su(ii),...
%        sqrt(abs(vu(ii))),sqrt(abs(su(ii))) ); 
    fprintf(fid,'%16.8e%16.8e%16.8e%16.8e\n', lon(jj),lat(jj),vu(jj),su(jj) ); 
end
fclose(fid);

% additional file sorted by decreasing MAGNITUDE

% % plot the scale to file
% filename = [filetag '_psxy_vert_scale.dat'];
% fid = fopen(filename,'w');
% fprintf(fid,'%18.6e', fscale);   
% fclose(fid);
    
%=======================================================================
