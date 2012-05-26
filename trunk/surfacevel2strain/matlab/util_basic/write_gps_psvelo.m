function write_gps_psvelo(filetag,lon,lat,ve,vn,se,sn,ren,name)
%WRITE_GPS_PSVELO writes a GPS vector field into psvelo format for GMT
%
% -Sevelscale/confidence/fontsize.
%       Velocity  ellipses  in  (N,E)  convention.   Vscale sets the scaling of the velocity arrows.  This scaling gives inches
%       (unless c, i, m, or p is appended).  Confidence sets the 2-dimensional confidence limit for the ellipse, e.g., 0.95 for
%       95%  confidence  ellipse.   Fontsize sets the size of the text in points.  The ellipse will be filled with the color or
%       shade specified by the -G option [default transparent].  The arrow and the circumference of the ellipse will  be  drawn
%       with the pen attributes specified by the -W option.  Parameters are expected to be in the following columns:
% 
% 1,2    longitude, latitude of station (-: option interchanges order)
% 3,4    eastward, northward velocity (-: option interchanges order)
% 5,6    uncertainty of eastward, northward velocities (1-sigma) (-: option interchanges order)
% 7      correlation between eastward and northward components
% 8      name of station (optional).
%

% number of stations
nstation = length(lon);

if nargin==8
   name = repmat(cellstr(''),nstation,1);
   for ii=1:nstation, name{ii} = sprintf('%4.4i',ii); end
end

if length(unique([length(lon) length(lat) length(ve) length(vn) ...
        length(se) length(sn) length(ren) length(name)]))~=1
   whos lon lat ve vn se sn ren name
   error('all input must have the same length');
end

filename = [filetag '_psvelo.dat'];
disp(['write_gps_psvelo.m: writing ' filename]);
fid = fopen(filename,'w');
for ii = 1:nstation
    fprintf(fid,'%12.4f%12.4f%12.4e%12.4e%12.4e%12.4e%12.4e%12s\n',...
        lon(ii),lat(ii),ve(ii),vn(ii),se(ii),sn(ii),ren(ii),char(name(ii)));   
end
fclose(fid);

% also write a file for psxy format (-Sv)
[th,r] = cart2pol(ve,vn);
filename = [filetag '_psxy.dat'];
disp(['write_gps_psvelo.m: writing ' filename]);
fid = fopen(filename,'w');
for ii = 1:nstation
    fprintf(fid,'%12.4f%12.4f%12.4e%12.4e%12s\n',...
        lon(ii),lat(ii),th(ii),r(ii),char(name(ii)));   
end
fclose(fid);    

%=======================================================================