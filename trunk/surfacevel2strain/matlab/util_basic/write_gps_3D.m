%
% function write_gps_3D.m
%
% This program writes a 3D GPS velocity field into a semi-standard format
% that can easily be loaded back into Matlab for other purposes.
%
% See also read_gps_3D.m for reading such files.
%
% calls xxx
% called by xxx
%

function write_gps_3D(filetag,lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name)

% number of stations
nstation = length(name);

filename = [filetag '_3D.dat'];
fid = fopen(filename,'w');
fprintf(fid,['%10s%10s' repmat('%12s',1,9) '%10s%10s%10s\n'],...
    'lon','lat','ve','vn','vu','se','sn','su','ren','reu','rnu','start','finish','name');
for ii = 1:nstation
    fprintf(fid,['%10.4f%10.4f' repmat('%12.4e',1,9) '%10.1f%10.1f%10s\n'],...
        lon(ii),lat(ii),ve(ii),vn(ii),vu(ii),...
        se(ii),sn(ii),su(ii),ren(ii),reu(ii),rnu(ii),...
        start_date(ii),finish_date(ii),name{ii} );   
end
fclose(fid);
    
%=======================================================================