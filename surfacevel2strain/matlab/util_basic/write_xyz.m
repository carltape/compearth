function write_xyz(filetag,x,y,z,stfmt)
%WRITE_XYZ write out three columns of x-y-z data

x = x(:); y = y(:); z = z(:);

filename = [filetag '_xyz.dat'];
n = length(x);

disp(sprintf(['write_xyz.m: ' filetag]));

if ~exist('stfmt','var')
    stfmt = '%18.10e%18.10e%18.10e';
end

% COLUMNS : x y z
fid = fopen(filename,'w');
for ii = 1:n
    fprintf(fid,[stfmt '\n'],x(ii),y(ii),z(ii) );   
end
fclose(fid);