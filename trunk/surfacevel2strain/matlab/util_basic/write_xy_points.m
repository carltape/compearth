%
% function write_xy_points(filetag,x,y)
%
% This file writes a simple 2-column file of x-y (or lon-lat) points.
%
% calls xxx
% called by xxx
%

function write_xy_points(filetag,x,y)

% number of points
n = length(x);

filename = [filetag '_points.dat'];
fid = fopen(filename,'w');
for ii = 1:n
    fprintf(fid,'%18.10e%18.10e\n',x(ii),y(ii) );   
end
fclose(fid);
    
%=======================================================================