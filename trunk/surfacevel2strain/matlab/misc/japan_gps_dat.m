%
% japan_gps_dat.m
% Carl Tape, 20-Mar-2009
%
% This file reads in the japan GPS dataset and makes some crude plots.
% The output files are written in a format that is ready to load for
% surfacevel2strain.m.
%
% calls read_gps_japan.m
% called by xxx
%

clear
close all
clc
format short, format compact

earthr = 6371*1e3;
earths = 4*pi*earthr^2;
deg = 180/pi;

%----------------------

iwrite = 0;   % =1 to write files

ofile = 'japan_takeo_ito';
odir = '/home/carltape/gmt/gps_data/ASIA/japan/';

% region of interest
isubregion = 1;             % =1 if you want to exclude observations outside this
ax0 = [128 147 30 46];

%=======================================

% read in the japan gps file
[lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] = read_gps_japan;
nobs = length(name);

% write the japan gps file for plotting in GMT
if iwrite == 1
    filetag = [odir ofile];
    
    % lon-lat points
    write_xy_points(filetag,lon,lat);
    
    % horizontal field
    write_gps_psvelo(filetag,lon,lat,ve,vn,se,sn,ren,name);
    
    % vertical field
    write_gps_psxy_vert(filetag,lon,lat,vu,su);

    % write for loading into Matlab
    write_gps_3D(filetag,lon,lat,ve,vn,vu,...
        se,sn,su,ren,reu,rnu,start_date,finish_date,name);
end

%=======================================

figure; hold on; plot(lon,lat,'.');

figure; hold on; quiver(lon,lat,ve,vn,'b');
plot(ax0([1 2 2 1 1]),ax0([3 3 4 4 3]),'r','linewidth',2);
axis equal
orient landscape, wysiwyg

% norm of the errors
err_norm = zeros(nobs,1);
for ii=1:nobs
    err_norm(ii) = norm([ se(ii) sn(ii) su(ii)]);
end

%------------------------------------

% fill data array
%data = [lon lat ve_res vn_res vu se sn su err_norm start_date];
data = [lon lat ve vn vu se sn su ren reu rnu err_norm start_date finish_date];

%------------------------------------
% ELIMINATE PARTICULAR OBSERVATIONS

if 0==1
    data5 = data;
    nobs5 = nobs;
    name5 = name;

else

    % check which velocities have sigma = 0 associated with one of the components
    izeroerr_e = find( se == 0); length(izeroerr_e)
    izeroerr_n = find( sn == 0); length(izeroerr_n)
    izeroerr_u = find( su == 0); length(izeroerr_u)

    izeroerr = unique([izeroerr_e ; izeroerr_n ; izeroerr_u ]);
    inonzero = setdiff([1:nobs]',izeroerr(:));
    %inonzero = vec_anti(nobs, izeroerr);

    plot(lon(izeroerr_e),lat(izeroerr_e),'ro');
    plot(lon(izeroerr_n),lat(izeroerr_n),'ro');
    plot(lon(izeroerr_u),lat(izeroerr_u),'ro');

    %figure; hold on; plot(ve,'b.'); plot(ve+se,'r.'); plot(ve-se,'r.'); orient landscape, wysiwyg
    %figure; hold on; plot(vn,'b.'); plot(vn+sn,'r.'); plot(vn-sn,'r.'); orient landscape, wysiwyg
    %figure; hold on; plot(vu,'b.'); plot(vu+su,'r.'); plot(vu-su,'r.'); orient landscape, wysiwyg

    % eliminate observations with sigma = 0
    data2 = data(inonzero,:);
    name2 = name(inonzero);
    nobs2 = length(name2);

    % find unique observation points
    [lonlat_unique,iind,jind] = unique(data2(:,1:2),'rows');
    nunique = length(iind);

    % loop over the unique points
    ikeep = [];
    for ii=1:nunique
        lon_temp = lonlat_unique(ii,1);
        lat_temp = lonlat_unique(ii,2);
        inds = find( and( data2(:,1) == lon_temp, data2(:,2) == lat_temp) );

        % sort by increasing norm of error, then take the first in the list
        [dtemp,dinds] = sortrows([ data2(inds,1) data2(inds,2) data2(inds,9) ],3);
        ikeep = [ikeep inds(dinds(1))];
        if length(inds) > 1
            name2(inds)
            disp(dtemp);
        end
    end
    data3 = data2(ikeep,:);
    name3 = name2(ikeep);
    nobs3 = length(name3);

    % extract observations within region of interest
    if isubregion==1
        inds = getsubset(data3(:,1),data3(:,2),ax0);
    else
        inds = [1:nobs3]';
    end
    data4 = data3(inds,:);
    name4 = name3(inds);
    nobs4 = length(name4);

    % sort from south to north
    [data5,i5] = sortrows(data4,2);
    name5 = name4(i5);
    nobs5 = length(name5);

    disp('  ');
    disp(sprintf('%30s%8i',' Initial data set: ',nobs));
    disp(sprintf('%30s%8i','Eliminating sigma = 0 obs: ',nobs2));
    disp(sprintf('%30s%8i','Eliminating redundant obs : ',nobs3));
    disp(sprintf('%30s%8i','Keeping obs within region : ',nobs4));
end

%------------------------------------

% write the japan gps field for plotting in GMT
% data = [lon lat ve vn vu se sn su ren reu rnu err_norm start_date finish_date];
if iwrite == 1
    % file label
    filetag = [odir ofile];
    
    % FULL japan field
    write_xy_points(filetag,lon,lat);
    write_gps_psvelo(filetag,lon,lat,ve,vn,se,sn,ren,name);
    write_gps_psxy_vert(filetag,lon,lat,vu,su);
    write_gps_3D(filetag,lon,lat,ve,vn,vu,...
        se,sn,su,ren,reu,rnu,start_date,finish_date,name);
    
    % REDUCED japan field
    % data = [lon lat ve vn vu se sn su err_norm start_date];
    lon5 = data5(:,1); lat5 = data5(:,2);
    ve5 = data5(:,3); vn5 = data5(:,4); vu5 = data5(:,5);
    se5 = data5(:,6); sn5 = data5(:,7); su5 = data5(:,8);
    ren5 = data5(:,9); reu5 = data5(:,10); rnu5 = data5(:,11);
    start_date5 = data5(:,13); finish_date5 = data5(:,14);
    %ren5 = zeros(nobs5,1); reu5 = zeros(nobs5,1); rnu5 = zeros(nobs5,1); finish_date5 = zeros(nobs5,1);
    
    filetag2 = [filetag '_subset'];
    
    write_xy_points(filetag2,lon5,lat5);
    write_gps_psvelo(filetag2,lon5,lat5,ve5,vn5,se5,sn5,ren5,name5);
    write_gps_psxy_vert(filetag2,lon5,lat5,vu5,su5);
    write_gps_3D(filetag2,lon5,lat5,ve5,vn5,vu5,...
        se5,sn5,su5,ren5,reu5,rnu5,...
        start_date5,finish_date5,name5);
end

%============================================================
