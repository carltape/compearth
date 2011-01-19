%
% function write_euler_vector_psxy.m
% CARL TAPE, 12-Aug-2007
% printed xxx
%
% This file writes a set of euler vectors to a file for plotting in GMT
% using psxy.
%
% Copied from write_gps_psxy_vert.m on 8-27-07.
%
% The sqrt is used so that the AREA of the plotting circle is proportional
% to the scalar value.
%
% The input scaling option controls the size of the markers.
%
% calls xxx
% called by xxx
%

function write_euler_vector_psxy(filetag,lon,lat,omega)

lon = lon(:);
lat = lat(:);
omega = omega(:);

% sort by DECREASING magnitude, so that the smaller dots are not covered
dall = [lon lat omega];
[dsort,isort] = sortrows(dall,[-3]);
lon = dsort(:,1);
lat = dsort(:,2);
omega = dsort(:,3);

nstation = length(lon);

% euler vectors (omega > 0)
filename = [filetag '_euler_vector_psxy.dat'];
fid = fopen(filename,'w');
for ii = 1:nstation
    %fprintf(fid,'%16.8e%16.8e%16.8e\n',lon(ii),lat(ii),omega(ii),abs(omega(ii))) );   
    fprintf(fid,'%16.8e%16.8e%16.8e\n',lon(ii),lat(ii),omega(ii) ); 
end
fclose(fid);

% anti-euler vectors (omega < 0)
[lat_anti,lon_anti] = antipode(lat,lon);   % requires Matlab mapping toolbox
filename = [filetag '_euler_anti_vector_psxy.dat'];
fid = fopen(filename,'w');
for ii = 1:nstation
    fprintf(fid,'%16.8e%16.8e%16.8e\n',lon_anti(ii),lat_anti(ii),-omega(ii) );   
end
fclose(fid);

% % plot the scale to file
% filename = [filetag '_euler_vector_scale_psxy.dat'];
% fid = fopen(filename,'w');
% fprintf(fid,'%18.6e', fscale);   
% fclose(fid);
 
%=======================================================================
