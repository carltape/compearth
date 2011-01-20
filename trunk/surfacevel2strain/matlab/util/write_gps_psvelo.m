%
% function write_gps_psvelo(filetag,lon,lat,ve,vn,se,sn,ren,name)
%
% This file writes a GPS vector field into a format compatible with the GMT
% plotting command psvelo.
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
% calls xxx
% called by xxx
%

function write_gps_psvelo(filetag,lon,lat,ve,vn,se,sn,ren,name)

% number of stations
nstation = length(name);

filename = [filetag '_psvelo.dat'];
fid = fopen(filename,'w');
for ii = 1:nstation
    fprintf(fid,'%12.4f%12.4f%12.4e%12.4e%12.4e%12.4e%12.4e%12s\n',...
        lon(ii),lat(ii),ve(ii),vn(ii),se(ii),sn(ii),ren(ii),char(name(ii)));   
end
fclose(fid);
    
%=======================================================================