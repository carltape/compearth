function [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] ...
    = platemodel2gps(lon,lat,imodel,ifix,opts)
%PLATEMODEL2GPS compute surface velocities from a set of plates and euler poles
%
% INPUT:
%   lon,lat     gridpoints for which you want the surface velocities
%   imodel      index into plate model (euler vectors and plate boundaries)
%   ifix        index of fixed plate (99 to leave vfield in the original ref frame)
%   opts        misc options
%
% OUTPUT:
%   lon,lat     gridpoints for which you have the surface velocities
%   ve, vn      surface velocities (mm/yr or km/Myr)
%   iplate_vec  indexing vector for the gridpoint
%                 iplate_vec(4) = 8 --> the 4th gridpoint is on plate number 8
%   exyz        euler vectors (wx,wy,wz), in deg/Myr; exyz(:,ifix) = [0 0 0]'
%   names       names of plates corresponding to the plate boundary files
%   name_labs   abbreviated labels of the plates
%
% imodel
%   (1)         Craig O'Neill model (from Eh Tan)
%   (2)         NUVEL-1A-nnr
%   (3)         REVEL
%   (4)         Bird
%   (5)         Gripp and Gordon (2002) -- HS3-NUVEL1A
%   (6)         Bird, but in GrippGordon reference frame
%   (7)         Bird, but in MorganMorgan reference frame
%   (8)         Bird, but in no-net-rotation reference frame
%
% The program will output the surface velocity vector field,
% as long as there are <5000 gridpoints.
%
% Originally, I adapted a pre-computed velocity field from Eh Tan -- those
% codes are in plate_model.m.  The key problem is to determine whether a
% particular point on the sphere is inside or outside of a plate that is
% described in terms of a closed-contour spherical polygon cap.
%
%==============================
% See notes in /home/carltape/latex/notes/misc/plate_models.pdf
%==============================
%
% EXAMPLES:
%    [lon,lat] = gridvec(-120,-114,50,32,37); platemodel2gps(lon,lat,2,11,{1,1,0});
%    [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] = platemodel2gps([],[],2,11,{0,1,1});
%
% calls:
%    get_plate_model.m
%    euler2gps.m
%    polygon_centroid_3d.m
%    euler_rot_tec.m
%    euler_convert.m
%
%    arcdist.m
%    unit.m
%    griddataXB.m
%    xyz2latlon.m, latlon2xyz.m
%
% called by run_platemodel2gps.m
%
% Carl Tape, 2006-01-11
%


% constants
deg = 180/pi;
earthr = 6371*1e3;

% check if input longitudes are [0,360] or [-180,180]
if any(lon > 180), ilon360 = 1; else ilon360 = 0; end

% might need to make sure that longitude is [-180,180]
lat = lat(:);
lon = lon(:);
num = length(lat);

% options
ifig_extra  = opts{1};
idisplay    = opts{2};
ieuler_only = opts{3};

% KEY: get euler vectors for plate model
%get_plate_model;   % input imodel
[exyz,names,name_labs,dir_bounds,ssfx,smod] = get_plate_model(imodel);
nump = length(names);

% display info on euler poles
if and(idisplay==1, ieuler_only==1)
    
    % convert euler poles (wx,wy,wz) --> (lat,lon,omg)
    outvec1 = euler_convert(exyz,1);
    elat = outvec1(1,:);
    elon = outvec1(2,:);

    % rotation rate, deg/Myr
    omegs = sqrt( exyz(1,:).^2 + exyz(2,:).^2 + exyz(3,:).^2 );

    % maximum surface velocity on a plate, in mm/yr
    % vmax will be at the points that are del=90 from the euler pole
    % vmax = omegs/mmyr2degMyr;
    
    % display info on euler poles
    disp('  '); disp(['Here are ' num2str(nump) ' euler poles from the ' smod ' model:']);
    disp('             label     elon        elat    deg/Myr    name'); 
    for ii=1:nump  
        disp(['  Plate #' num2str(sprintf('%2.2i',ii)) ' : ' ...
            name_labs{ii} ...
            num2str(sprintf('%12.4f',elon(ii))) ...
            num2str(sprintf('%12.4f',elat(ii))) ...
            num2str(sprintf('%10.4f',omegs(ii))) ...
            '   ' names{ii} ])
    end
    disp('-----------------------------');
end

disp('  ');
if ifix==99
    disp(' will not change reference field for the velocity field');
    flab = ' in unchanged reference frame';
else
    % fix one of the plates by subtracting its euler vector from all vectors
    disp([' fixed plate is ' name_labs{ifix} ' (' names{ifix} ')']);
    flab = [' with fixed ' name_labs{ifix}];
    exyz0 = exyz;
    exyz = exyz0 - repmat( exyz0(:, ifix), 1, nump );
end

% if all you want is the plate model, then exit
if ieuler_only==1
    ve = []; vn = []; iplate_vec = [];
    return
end

%========================================================
% ASSIGN PLATE INDEX TO EACH GRIDPOINT
%   In essence, this involves rotating each plate so that its centroid is at
%   (lat=0, lon=0). The same finite rotation is applied to the test
%   gridpoints. Then we use the Matlab command inpolygon.m to determine
%   whether a gridpoint is inside the boundary.
%
%   This is not an ideal algorithm, since inpolygon.m assumes edges that are
%   not arcs, but rather chords in the Cartesian lat-lon domain.  This
%   problem would be overcome with extremely dense sampling of the plate
%   boundaries.  But for regional plotting purposes, using inpolygon.m is
%   fine.

if length(lon)==1
    ax1 = [lon-1 lon+1 lat-1 lat+1];
else
    ax1 = [min(lon) max(lon) min(lat) max(lat)];
end

iplate_vec = zeros(num,1);

disp('platemodel2gps.m: find the plate for each target point');

% loop over plates
pmin = 1; pmax = nump;
for ii = pmin:pmax
    disp(sprintf('plate is %s, index %i (%i to %i)',names{ii},ii,pmin,pmax));
    
    % load plate boundary file (assumes longitude in column 1)
    if any(imodel==[9 10]), ww = name_labs{ii}; else ww = names{ii}; end
    load([dir_bounds ww ssfx]);
    data_plot = eval(ww);
    plon = data_plot(:,1);
    plat = data_plot(:,2);
    
    %-----------------------
    % determine finite rotation so that centroid is at (plat,plon) = (0,0)
    
    % xyz points of the plate boundary
    Pxyz = latlon2xyz(plat,plon,earthr);
    n = length(plat);
    
    % (approximate) centroid of the spherical patch
    centroid = polygon_centroid_3d(n, Pxyz);
    [clat,clon] = xyz2latlon(centroid);
    
    % angle of finite rotation
    Zxyz = [earthr 0 0]';
    [Zlat,Zlon] = xyz2latlon(Zxyz);
    dist = arcdist(clat, clon, Zlat, Zlon);
    
    % pole of finite rotation
    [Aunit, Alen] = unit(cross(centroid,Zxyz) );
    [Elat,Elon] = xyz2latlon(Aunit);
    evec = [Elat Elon dist]';
    
    % rotate plate boundary
    [plat_rot, plon_rot] = euler_rot_tec(plat,plon,evec);
    
    % rotate gridpoints
    [lat_rot, lon_rot] = euler_rot_tec(lat,lon,evec);
    
    %-----------------------
    % un-rotate the rotated plate to check (note -dist)
    Pxyz = latlon2xyz(plat_rot,plon_rot,earthr);
    [plat_orig, plon_orig] = euler_rot_tec(plat_rot, plon_rot, [Elat Elon -dist]');
    
    %-----------------------
    % determine which gridpoints are inside the plate
    % (we assume none of the gridpoints are EXACTLY on the boundaries)
    
    xv = plon_rot; yv = plat_rot;
    x = lon_rot; y = lat_rot;
    in = inpolygon(x,y,xv,yv);          % KEY COMMAND
    
    iplate_vec(in) = ii;
    
    if ifig_extra==1
        figure; hold on;
        plot(lon(in),lat(in),'.k');
        plot(plon,      plat,'r');
        plot(plon_rot,  plat_rot,'b');
        plot(plon_orig, plat_orig,'k--');
        plot(clon,clat,'r.','markersize',16);
        plot(Elon,Elat,'g.','markersize',16);
        plot(Zlon,Zlat,'b.','markersize',16);
        axis equal, axis([-180 180 -90 90]); title([num2str(ii) ' : ' ww]);
        fontsize(9), orient tall, wysiwyg
    end
end

if ifig_extra==1
    figure; hold on;
    [X,Y,Z] = griddataXB(lon,lat,iplate_vec,400,'nearest');
    pcolor(X,Y,Z); shading flat; caxis([pmin pmax])
    %plot(lon, lat,'k+');
    axis equal, axis(ax1);
    colorbar

    pmin = 1; pmax = nump;
    for ii = pmin:pmax   

        % CHECK: load data and plot
        ww = names{ii};
        load([dir_bounds ww ssfx]);
        data_plot = eval(ww);
        plon = data_plot(:,1);
        plat = data_plot(:,2);

        plot(plon, plat,'k.');
    end
end

%========================================================
% COMPUTE SURFACE VELOCITY FIELD

disp('platemodel2gps.m: compute surface velocity fields');

vn = zeros(num,1);
ve = zeros(num,1);
for ivel = 1:nump

    % gridpoints on plate ivel
    inds = find( ivel==iplate_vec );

    % compute surface velocity field
    evec = exyz(:,ivel);                            % euler pole
    Pxyz = latlon2xyz(lat(inds),lon(inds),earthr);	% input points
    Vrtp = euler2gps(evec, Pxyz);                   % surface vel (local r,th,ph)
    vn(inds) = -Vrtp(2,:)';
    ve(inds) = Vrtp(3,:)';
end

% plot vector field (as long as there are not too many vectors)
if num <= 5000
    figure; hold on;
    quiver(lon,lat,ve,vn,1);
    title(sprintf('platemodel2gps.m: model %s %s (%i plates, %i gridpoints)',...
        smod,flab,nump,num),'interpreter','none');
    xlabel('Longitude'); ylabel('Latitude');
    axis equal, axis(ax1);
    
    pmin = 1; pmax = nump;

    % plot plate boundaries
    for ii = pmin:pmax   
        if any(imodel==[9 10]), ww = name_labs{ii}; else ww = names{ii}; end
        load([dir_bounds ww ssfx]);
        data_plot = eval(ww);
        plon = data_plot(:,1);
        if ilon360==1, plon = wrapTo360(plon); end
        plat = data_plot(:,2);

        plot(plon, plat,'k.');
    end
else
    disp(sprintf('platemodel2gps.m: %i points, so not plotting vector field',num));
end

disp('returning from platemodel2gps.m');

%========================================================
