function write_xy_points(filetag,x,y,stfmt)
%WRITE_XY_POINTS writes a simple 2-column file of x-y (or lon-lat) points

% number of points
n = length(x);

if ~exist('stfmt','var')
    stfmt = '%18.10e%18.10e';
end

filename = [filetag '_points.dat'];
disp(sprintf('writing %i points to file %s',n,filename));
fid = fopen(filename,'w');
for ii = 1:n
    fprintf(fid,[stfmt '\n'],x(ii),y(ii) );   
end
fclose(fid);
