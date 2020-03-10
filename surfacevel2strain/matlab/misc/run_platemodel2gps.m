%
% run_platemodel2gps.m
%
% This program computes the surface velocity field at a set of points on
% the sphere. The Euler vectors and plate boundaries are provided from
% plate models. A key problem is to determine whether a particular point
% on the sphere is inside or outside of a plate that is
% described in terms of a closed-contour spherical polygon cap.
%
% EXAMPLES (note: PA is fixed in the Bird plate model):
% (1) plate velocity vectors shown in Wang and Tape (2014), Figure 11a
%     DEFAULT OPTIONS        // SUBSEQUENT EXAMPLES
%     4 (plate model = bird) // 8, 7, 6
%     1 (grid)               
%     4 (q)                  
%     2 (npac)               
%     1 (write)               
%     1 (fixed plate)        // 0, 0, 0
%     11 (fixed NAM)         // 
% (2) user input for manually entered lon-lat points
%     4 (plate model = bird) // 6, 7
%     3 (manually entered points)                    
%     2 (npac)               
%     0 (no write)               
%     0 (no fixed plate)
%
% EXAMPLE 2 OUTPUT:
% surface velocities (lon, lat, ve (mm/yr), vn (mm/yr)):
%    -148.4460      70.2036     3.20   -15.29
%    -149.5900      68.6274     3.03   -15.11
%    -150.1100      67.3812     2.79   -15.03
%    -150.2640      66.2066     2.49   -15.01
%    -148.5110      65.5114     1.86   -15.28
%
% (3) figure shown in plate_model.pdf notes
%     7 (plate model = bird_morgan)
%     1 (grid)
%     7 (q)                   // 7
%     1 (globe)
%     1 (write) 
%     0 (fixed plate)
%
% To calculate surface velocities for points that are outside the plate,
% e.g., above subducting plates, use platemodel2conv_vel.m.
% 
% calls platemodel2gps.m
%
% Carl Tape, 2011-01-28
%

clear, clc, close all
format short
format compact

% add path to additional matlab scripts
user_path;
%========================================================
% USER PARAMETERS

% plate models
mod_labs = {'oneill','nuvel1A_nnr','revel','bird','gripp_hs3','bird_gripp','bird_morgan','bird_nnr','','morvel_nnr'};
nmod = length(mod_labs);
disp(sprintf('\nPLATE MODELS (INCLUDING REFERENCE FRAME) TO CHOOSE FROM:'));
for ii=1:nmod, disp(sprintf('%5i  %s',ii,mod_labs{ii})); end
disp('11-21  morvel_nnr with a reference frames in Becker et al. (2015)');
imodel = input(sprintf(' Type index for plate model (1-%i), then ENTER: ',nmod));

% type of grid points
disp('  '); disp('PLOTTING GRID POINTS TO CHOOSE FROM:');
disp('  igrid = 1 -- use spherical-triangle grids (q = 0,1,..)');
disp('  igrid = 2 -- use regular mesh of lon-lat points');
disp('  igrid = 3 -- use arbitrary points');
igrid = input(' Type index for grid points (1,2,3), then ENTER: ');

% parameters controlling the number of gridpoints
if igrid==1
    q = input(' Type spherical grid index (q=0,1,2,...), then ENTER: ');
elseif igrid==2
    q = 99;
    numx = input(' Type grid spacing along x (lon) direction, then ENTER: ');
elseif igrid==3
    q = 99;
    disp(' load lon-lat points next');
end

% user-specified regions
sregions = {'globe','npac','alaska','socal','japan','jdf','micros','sumatra',...
    'taiwan','tibet','eurasia','yakutat'};
nregion = length(sregions);
disp(sprintf('\nGEOGRAPHIC REGIONS TO CHOOSE FROM:'));
for ii=1:nregion
    disp(sprintf('%3i  %s',ii,sregions{ii}));
end
iregion = input(sprintf(' Type index for region (1-%i), then ENTER: ',nregion));
switch iregion
    case 1, ax1 = [-180 180 -90 90];
    case 2, ax1 = [154 250 30 74];
    case 3, ax1 = [154 250 30 74];
    case 4, ax1 = [-122 -113 30 38];
    case 5, ax1 = [100 160 5 60];
    case 6, ax1 = [-135 -115 35 55];
    case 7, ax1 = [115 160 -15 20];
    case 8, ax1 = [85 110 -10 30];
    case 9, ax1 = [115 125 18 28];
    case 10, ax1 = [52 104 12 44];
    case 11, ax1 = [-10 125 5 50];  
    case 12, ax1 = [-147 -136 57 63];    
end
slabel = sregions{iregion};
disp(sprintf('region %i is for %s: [%.1f %.1f %.1f %.1f]',iregion,slabel,ax1));

iwrite = input(' Type 1 to write files (0 otherwise), then ENTER: ');
ifig_extra = 0;     % extra figures
ipick_figs = 0;     % figures of example plate w.r.t. fixed plate
ilon360 = 1;        % =1 for longitudes as [0,360], =0 for [-180,180]

stq = sprintf('q%2.2i', q);
sgrid = sprintf('g%1i',igrid);
if imodel <= 10
    smod = mod_labs{imodel};
else
    smod = sprintf('morvel_nnr_becker2015_%2.2i',imodel-10);
end
dir_plates  = '/home/carltape/gmt/plates/';
dir_surface = [dir_plates 'surface_velocities/' smod '/'];
dir_models  = [dir_plates 'plate_models/' smod '/'];
%dir_grids   = '/home/carltape/SURFACEVEL2STRAIN/fortran/grids_output/full_grids/';

%==========================================================================

if ilon360==1
    lontick = [0:60:360];
else
    lontick = [-180:60:180];
end
lattick = [-90:30:90];
deg = 180/pi;
earthr = 6371*1e3;
msize = 10^2;

% convert mm on earth surface to deg
mm2deg = 1e-3 / 6371e3 * deg;

% convert mm/yr on earth surface to deg/Myr
mmyr2degMyr = mm2deg * 1e6;

%==========================================================================
% INPUT GRIDPOINTS FOR VELOCITY VECTORS AND BACKGROUND VELOCITY MAGNITUDE

lonmin = ax1(1); lonmax = ax1(2);
latmin = ax1(3); latmax = ax1(4);

if igrid==1

    [lon,lat,~,~,~] = getspheregrid(ax1,q,q);  
    [~,lon,lat] = getsubset(lon,lat,ax1);

elseif igrid == 2
    
    [lon,lat] = gridvec(lonmin,lonmax,numx,latmin,latmax);

%     % denser uniform mesh for plotting the magnitude of the velocity field
%     numx_mag = 10*numx;
%     xvec = linspace(lonmin,lonmax,numx_mag);
%     dx_mag = xvec(2) - xvec(1);
%     xvec = xvec(1:end-1); numx_mag = length(xvec);
%     dy_mag = dx_mag;
%     yvec = [latmin : dy_mag : latmax];
%     numy_mag = length(yvec);
%     [X,Y] = meshgrid(xvec,yvec);
%     lon_mag = X(:);
%     lat_mag = Y(:);

elseif igrid==3
    if 0==1
        % click points that you want
        [lon,lat] = getxy(ax1);
    else
        % manually list points (or load from file)
        disp('USING MANUALLY ENTERED POINTS');
        lon = [-148.446 -149.59 -150.11 -150.264 -148.511];
        lat = [70.2036 68.6274 67.3812 66.2066 65.5114];
    end
end

% exclude the poles, since you have no concept of east and north there
ipoles = find(abs(lat) > 89.9);
lon(ipoles) = [];
lat(ipoles) = [];

num = length(lat);

%==========================================================================
% EULER POLES --- ALL PLATE MODELS
% See notes in /home/carltape/latex/notes/misc/plate_models.pdf

if 0==1
    for imodel = 1:7
        platemodel2gps([],[],imodel,99,{0,1,1});
    end
    disp('displaying euler poles only');
    error
end

%==========================================================================
% EULER POLES FOR PLATE MODEL (see test_euler2gps.m)
% THESE ARE PRE-SAVED AND THEN LOADED BY get_plate_model.m FROM INSIDE platemodel2gps.m

dir_models = [dir_plates 'plate_models/' smod '/'];
ww = [dir_models smod '_euler_poles.dat'];

imake = 0;  % construct the (wx,wy,wz) data file for the plate model

% for imodel = 11:21
%    smod = sprintf('morvel_nnr_becker2015_%2.2i',imodel-10);
%    dir_models  = [dir_plates 'plate_models/' smod '/'];
%    ww = [dir_models smod '_euler_poles.dat'];
%    imodel
    
if imake==1
    % load ORIGINAL euler poles and output in common format: wx, wy, wz, LAB
    switch imodel
        case 1
            % load rotation vectors (units in deg/Myr)
            % (see test_euler_rot_tec.m)
            [name_labs,wx,wy,wz,names] = textread([dir_models 'ind2002_mod.a1'],'%s%f%f%f%s','headerlines',2);
            exyz = [wx wy wz]';
            
            % indexing into euler poles used in making the plate velocity model
            % (take 13 out of 22 plates)
            inds = [3 1 10 2 9 8 5 14 7 6 4 11 12];
            exyz = exyz(:,inds);

            % three-letter names
            names = names(inds);
            name_labs = name_labs(inds);
            %name_labs = {'AFR','ANT','ARA','AUS','CAR','COC','EUR','JDF','NAZ','NAM','PAC','PHI','SAM'}';
            
        case 2
            % load rotation vectors (units in deg/Myr)
            % (see test_euler_rot_tec.m)
            [name_labs,elon,elat,omeg,wx,wy,wz,names] = textread([dir_models 'nnr_nuvel1a_mod.dat'],'%s%f%f%f%f%f%f%s','headerlines',1);
            exyz = [wx wy wz]';
            
            % convert rad/Myr --> deg/Myr
            [Eunit, Elen] = unit(exyz);
            exyz = Eunit .* repmat( Elen*deg, 3, 1 );
            
            % indexing into euler poles used in making the plate velocity model
            inds = [1:8 13 10 9 11 14 16 12];
            exyz = exyz(:,inds);
            names = names(inds);
            
            % three-letter names
            name_labs = {'AFR','ANT','ARA','AUS','CAR','COC','EUR','IND','JDF','NAZ','NAM','PAC','PHI','SCO','SAM'}';
            
        case 3
            % load rotation vectors (units in deg/Myr)
            % (see test_euler_rot_tec.m)
            % MODIFIED TO MATCH THE BIRD LABELING SCHEME
            [name_labs,wx,wy,wz,names] = textread([dir_models 'revel_mod.dat'],'%s%f%f%f%s','headerlines',1);
            exyz = [wx wy wz]';
            
            % convert 10^-3 rad/Myr --> deg/Myr
            [Eunit, Elen] = unit(exyz);
            exyz = Eunit .* repmat( Elen*10^-3*deg, 3, 1 );
            
            % indexing into euler poles used in making the plate velocity model
            %inds = [11 2 3 5 6 12 8 10 10 12 14 15 16];
            inds = [11 1 4 2 3 5 6 8 9 12 10 13 14 15 17 16 19 7];
            
            exyz = exyz(:,inds);

            % two-letter names (Bird plates)
            name_labs = name_labs(inds);
            names = names(inds);
            %name_labs = {'AF','AM','AT','AN','AR','AU','CA','EU','IN','NZ','NA','OK','PA','PS','SO','SA','SU','YA'}';
            %names = {'africa','amur','anatolia','antarctica','arabia','australia',...
            %    'caribbean','eurasia','india','nazca','north_america','okhotsk','pacific',...
            %    'philippine_sea','somalia','south_america','sunda','yangtze'}';
            
         case 4
            % load rotation vectors (units in deg/Myr)
            % (see test_euler_rot_tec.m)
            [name_labs,elat,elon,omega,names] = textread([dir_models 'PB2002_poles_mod2.dat'],'%s%f%f%f%s','headerlines',1);
            elatlon = [elat elon omega]';
            
            % convert to wx,wy,wz
            exyz = euler_convert(elatlon, 0);
            
            % indexing into euler poles used in making the plate velocity model
            %inds = [1 3 5 8 13 15 18 22 29 34 37 39 41];
            inds = 1:length(names);
            exyz = exyz(:,inds);
            names = names(inds);

            % three-letter names
            %name_labs = {'AFR','ANT','ARA','AUS','CAR','COC','EUR','JDF','NAM','NAZ','PAC','PHI','SAM'}';   
            
        case 5
            % load euler vectors (units in deg/Myr)
            [names,elat,elon,omega,name_labs] = textread([dir_models 'gripp_gordon_mod.dat'],'%s%f%f%f%s','headerlines',1);
            elatlon = [elat elon omega]';
            
            % convert to wx,wy,wz
            exyz = euler_convert(elatlon, 0);
            
            % indexing into euler poles used in making the plate velocity model
            %inds = [1 3 5 8 13 15 18 22 29 34 37 39 41];
            inds = 1:length(names);
            exyz = exyz(:,inds);
            names = names(inds);
            
        case 6
            % Bird model in hotspot reference frame, obtained by using the
            % NUVEL1A-HS3 Pacific plate euler vector
            
            % load the NUVEL1A-HS3 euler vectors
            [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] = platemodel2gps([],[],5,99,{0,1,1});
            exyz_pac_hs3 = exyz(:,12);
            
            % load the Bird euler vectors
            [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] = platemodel2gps([],[],4,99,{0,1,1});
            exyz_pac_bird = exyz(:,37);
            
            % Bird euler vectors in NUVEL1A-HS3 reference frame
            exyz = exyz - repmat(exyz_pac_bird - exyz_pac_hs3, 1, length(names));

            % indexing into euler poles used in making the plate velocity model
            inds = 1:length(names);
            exyz = exyz(:,inds);
            names = names(inds);
            
        case 7
            % Bird model in hotspot reference frame, obtained by using the
            % Morgan and Morgan (2007) Pacific plate euler vector
            
            % Morgan and Morgan (2007) Pacific plate euler vector
            % WHY DO I NEED TO FLIP THE SIGN ON THEIR LOCATION?
            exyz_pac_morgan = euler_convert([-59.44 84.91 0.8029],0);
            
            % load the Bird euler vectors
            [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] = platemodel2gps([],[],4,99,{0,1,1});
            exyz_pac_bird = exyz(:,37);
            
            % Bird euler vectors in NUVEL1A-HS3 reference frame
            exyz = exyz - repmat(exyz_pac_bird - exyz_pac_morgan, 1, length(names));
            
            % indexing into euler poles used in making the plate velocity model
            inds = 1:length(names);
            exyz = exyz(:,inds);
            names = names(inds);
            
        case 8
            % Bird model in no-net-rotation frame, obtained by using the
            % Morgan and Morgan (2007) Pacific plate euler vector
            
            % load the NNR-NUVAL1A euler vectors
            imodel0 = 2;
            [lon0,lat0,ve0,vn0,iplate_vec0,exyz0,names0,name_labs0] = platemodel2gps([],[],imodel0,99,{0,1,1});
            ipac0 = 12;
            exyz_pac_nnr = exyz0(:,ipac0);
            
            % load the Bird euler vectors
            [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] = platemodel2gps([],[],4,99,{0,1,1});
            exyz_pac_bird = exyz(:,37);
            
            % Bird euler vectors in NNR reference frame
            % note: subtract the Bird-Pacific, then add the NNR-Pacific
            exyz = exyz - repmat(exyz_pac_bird - exyz_pac_nnr, 1, length(names));
            
            % indexing into euler poles used in making the plate velocity model
            inds = 1:length(names);
            exyz = exyz(:,inds);
            names = names(inds);
            
        case 10
            [name_labs,elat,elon,omega,names] = textread([dir_models 'morvel_nnr.txt'],'%s%f%f%f%s');
            elatlon = [elat elon omega]';
            exyz = euler_convert(elatlon, 0);
            inds = 1:length(names);
            
        otherwise
            [lon,lat,ve,vn,iplate_vec,exyz,names,name_labs] = platemodel2gps([],[],10,99,{0,1,1});
            disp('using the NNR-MORVEL56 plate model');
            %iepick = input(sprintf(' Type index for reference frame from Becker2015 [0 = default morvel_nnr] [0-11]: ',nmod));
            iepick = imodel - 10;
            if iepick > 0
                [ename,elon,elat,omega] = textread('/home/carltape/papers/plate_models/Becker2015/Becker2015_Table1.txt','%s%f%f%f');
                elatlon = [elat elon omega]';
                exyz0 = euler_convert(elatlon,0);
                exyz = exyz + repmat(exyz0(:,iepick),1,length(names));
            end
            inds = 1:length(names);
            
    end

    % check the indexing
    disp([names name_labs])
    
    % write euler poles to text file
    % NEED TO MAKE SURE THAT THE DIRECTORY IS THERE!
    fid = fopen(ww,'w');
    for ii=1:length(inds)
        fprintf(fid,'%13.6f%13.6f%13.6f%16s%16s\n',exyz(1,ii),exyz(2,ii),exyz(3,ii),name_labs{ii},names{ii});   
    end
    fclose(fid);
    
    error
end   % imake==1

%end  % for imodel=11:21
%error

%==========================================================================
% COMPUTE SURFACE VELOCITY FIELD

disp('computing the surface velocity field...');

% INDEX OF THE FIXED PLATE IN EACH PLATE MODEL
% 99 indicates that the plate does not have an euler vector in the model
sfix_name = {'AFR','ANT','ARA','AUS','CAR','COC','EUR','IND','JDF','NAZ','NAM','PAC','PHI','SAM','SCO','SU','YA'};
ifix_mat = [ 1 1 1 2 1        %  1 AFR
             2 2 4 6 2        %  2 ANT
             3 3 5 7 3        %  3 ARA
             4 4 6 8 4        %  4 AUS
             5 5 7 13 5       %  5 CAR
             6 6 99 15 6      %  6 COC
             7 7 8 18 7       %  7 EUR
             99 8 9 21 8      %  8 IND
             8 9 99 22 9      %  9 JDF
             9 10 11 32 10    % 10 NAZ
             10 11 11 32 11   % 11 NAM
             11 12 13 37 12   % 12 PAC
             12 13 14 39 13   % 13 PHI (PS)
             13 15 16 46 14   % 14 SAM
             99 14 99 42 15   % 15 SCO
             99 99 17 48 99   % SU
             99 99 18 52 99   % YA
             ];
% the Bird plates are used in plate models 6, 7, 8
ifix_mat(:,6) = ifix_mat(:,4);
ifix_mat(:,7) = ifix_mat(:,4);
ifix_mat(:,8) = ifix_mat(:,4);
% the MORVEL56 plates are used in plate models 10,11-21
ifix_morvel = [14 3 4 5 6 7 9 10 11 16 13 1 17 19 20 23 25]';
% number of options for choosing the fixed plate
nfix = length(ifix_mat);

disp(' ');
disp('Type 0 to NOT fix a plate, 1 to fix a plate,');
ifix0 = input('or 2 to consider a range of fixed plates, then ENTER: ');
if ifix0 > 0
   disp('Index for fixed plate:');
   for ii=1:nfix, disp(sprintf('%3i %s',ii,sfix_name{ii})); end
   if ifix0 == 1
       ifix0 = input(sprintf('\n Type index of plate (between 1 and %i), then ENTER: ',nfix));
       irow1 = ifix0;
       irow2 = ifix0;
   else
       irow1 = input(sprintf('\n Type index of first fixed plate (between 1 and %i), then ENTER: ',nfix));
       irow2 = input(sprintf('\n Type index of first fixed plate (between %i and %i), then ENTER: ',irow1,nfix));
   end
   
else
    disp('Fixed reference plate is NOT being used');
    irow1 = 11; irow2 = irow1;
end

for irow = irow1:irow2     % KEY: loop over fixed plates
    
    if ifix0==0
        ifix = 99;
    else
        if imodel >= 10
            ifix = ifix_morvel(irow);
        else
            ifix_vec = ifix_mat(irow,:);
            ifix = ifix_vec(imodel);
        end
    end
    
    opts = {0,1,0};

    % compute velocity field for VECTORS (coarse mesh)
    % note: plots a vector field using the quiver command
    [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] ...
        = platemodel2gps(lon,lat,imodel,ifix,opts);
    
    if ifix==99, stref = '99'; else stref = name_labs{ifix}; end
    nump = length(names);

    % magnitudes: sort to plot largest magnitudes last (in GMT)
    vmag = sqrt(ve.^2 + vn.^2);
    
    % plot magnitudes as colored circles
    figure; hold on; grid on;
    scatter(lon,lat,msize,vmag,'filled');
    %scatter(lon,lat,msize,'ko');
    axis(ax1); colorbar;
    
    %======================================================================
    % WRITE VELOCITY FIELD TO FILE

    if iwrite==0
        disp('surface velocities (lon, lat, ve (mm/yr), vn (mm/yr)):');
        for ii=1:num
            disp(sprintf('%12.4f %12.4f %8.2f %8.2f',lon(ii),lat(ii),ve(ii),vn(ii)));
        end
        
    else
        % kludge for GMT plotting: set near-zero values to <0
        veps = 1e-4; vmag(vmag < veps) = -veps;
        % do not plot points that are fixed
        izero = find(vmag < veps);
        warning(sprintf('removing %i/%i points with v < %.3e mm/yr',length(izero),num,veps));
        lon(izero) = []; lat(izero) = []; ve(izero) = []; vn(izero) = []; vmag(izero) = [];
        num = length(lon);
        
        % output directory
        odir = [dir_plates 'surface_velocities/'];
        disp('writing to output directory: ');
        disp(odir);
        disp('columns are: lon lat ve vn vmag');
        
        % KEY: tag for all files
        ftag0 = [slabel '_fix_' stref '_' smod];
        ftag = [ftag0 '_' sgrid stq];
        
        % write plate model vector COMPONENTS and MAGNITUDES to file (mm/yr)
        [~,isort] = sort(vmag);
        ww = [ftag '_vec.dat'];
        ofile = [odir ww];
        disp(['writing file ' ww]);
        fid = fopen(ofile,'w');
        for ii=1:num
            jj = isort(ii);
            fprintf(fid,'%12.4f%12.4f%12.4f%12.4f%12.4f\n',...
                lon(jj),lat(jj),ve(jj),vn(jj),vmag(jj));   
            %fprintf(fid,'%18.8e%18.8e%18.8e%18.8e%18.8e\n',...
            %    lon(jj),lat(jj),ve(jj),vn(jj),vmag(jj));   
        end
        fclose(fid);
        
        % check
        if 0==1
            a = load('/home/carltape/gmt/plates/surface_velocities/globe_fix_99_bird_morgan_g1q03_vec.dat');
            figure; quiver(wrapTo360(a(:,1)),a(:,2),a(:,3),a(:,4));
        end

        % write bounds to file
        ww = [ftag '_bounds.dat'];
        ofile = [odir ww];
        disp(['writing to file ' ww]);
        fid = fopen(ofile,'w');
        fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',lonmin,lonmax,latmin,latmax); 
        
        %-------------------------------------
        % write euler poles
        % note: these files will be the same for many different cases
        
        % compute euler poles (lon-lat) and anti-poles (lon-lat)
        [elat,elon] = xyz2latlon(exyz);
        [elat_anti,elon_anti] = antipode(elat,elon);

        % write euler poles to file (lon-lat)
        ww = [ftag0 '_euler_lonlat.dat'];
        ofile = [odir ww];
        disp(['writing file ' ww]);
        fid = fopen(ofile,'w');
        for ii=1:nump, fprintf(fid,'%14.5f %14.5f \n',elon(ii),elat(ii) ); end
        fclose(fid);

        % write euler anti-poles to file (lon-lat)
        ww = [ftag0 '_euler_anti_lonlat.dat'];
        ofile = [odir ww];
        disp(['writing file ' ww]);
        fid = fopen(ofile,'w');
        for ii=1:nump, fprintf(fid,'%14.5f %14.5f \n',elon_anti(ii),elat_anti(ii) ); end
        fclose(fid);
        %-------------------------------------
        
    end

%==========================================================================
% CODE BELOW NEEDS TO BE REVISED

error

    % write to file for GMT plotting -- NON-GLOBAL FIELD (q = 99)
    if and(iwrite == 1, q == 99)

        % compute velocity field for MAGNITUDES (dense mesh)
        [lon_mag, lat_mag, ve_mag, vn_mag] = platemodel2gps(lon_mag,lat_mag,imodel,ifix,opts);
        vmag = sqrt( ve_mag.^2 + vn_mag.^2 );
        %plot(lon_mag,lat_mag,'r.');

        % for GMT plotting, set zero values to <0
        veps = 1e-4; vmag(vmag < veps) = -veps;

        % output directory
        odir = [dir_plates 'surface_velocities/misc/'];
        
        % write plate model vector COMPONENTS (fixed ifix plate) to file (mm/yr)
        ww1 = [slabel '_fix_' stref '_' smod '_vec.dat'];
        ofile = [odir ww1];
        disp(['writing to file ' ofile]);
        fid = fopen(ofile,'w');
        for ii=1:num
            fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',lon(ii),lat(ii),ve(ii),vn(ii));   
        end
        fclose(fid);

        % write plate model vector MAGNITUDES (fixed ifix plate) to file (mm/yr)
        ww2 = [slabel '_fix_' stref '_' smod '_mag.dat'];
        ofile = [odir ww2];
        disp(['writing to file ' ofile]);
        fid = fopen(ofile,'w');
        for ii=1:num_mag
            fprintf(fid,'%18.8e%18.8e%18.8e\n',lon_mag(ii),lat_mag(ii),vmag(ii));   
        end

        % write bounds to file
        ww3 = [slabel '_fix_' stref '_' smod '_bounds.dat'];
        ofile = [odir ww3];
        disp(['writing to file ' ofile]);
        fid = fopen(ofile,'w');
        fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',lonmin,lonmax,latmin,latmax); 
        if igrid==2
            fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',min(xvec),max(xvec),min(yvec),max(yvec)); 
            fprintf(fid,'%18.8e%18.8e\n',dx_mag,dy_mag); 
            fprintf(fid,'%14i%14i\n',numx_mag,numy_mag);
        end
    end
    
    %=======================================================
    % WRITE FIELDS TO FILE FOR GMT (/home/carltape/gmt/globe_figs/)

    % write the global velocity field and euler poles to file
    if and(iwrite == 1, q ~= 99)
        disp('  ');
        if ifix ~= 99
            disp([' fixed plate is ' names{ifix}]); 
        else
            disp(' no fixed plate: use original Euler vectors');
        end

        % for GMT plotting, set zero values to <0
        vmag = sqrt( ve.^2 + vn.^2 );
        veps = 1e-4; vmag(vmag < veps) = -veps;

        if q <= 4
            % write plate model vector COMPONENTS (fixed ifix plate) to file (mm/yr)
            ww = [smod '_fix_' stref '_vec_q' stq '.dat'];
            ofile = [dir_surface ww];
            disp(['writing to file ' ofile]);
            fid = fopen(ofile,'w');
            for ii=2:num-1  % exclude the poles
                fprintf(fid,'%14.5f %14.5f %12.4f %12.4f \n',lon(ii),lat(ii),ve(ii),vn(ii));   
            end
            fclose(fid);
        end

        % write plate model vector MAGNITUDES (fixed ifix plate) to file (mm/yr)
        vsort = sortrows([lon lat vmag],3);  % sort to plot largest magnitudes last
        ww = [smod '_fix_' stref '_mag_q' stq '.dat'];
        ofile = [dir_surface ww];
        disp(['writing to file ' ofile]);
        fid = fopen(ofile,'w');
        for ii=1:num
            fprintf(fid,'%14.5f%14.5f%12.4f \n',vsort(ii,:));   
        end
%         % repeat the lon=0 points at lon=360 for GMT plotting with grdimage
%         ilon0 = find(abs(lon) < 0.01);
%         for ii=1:length(ilon0)
%             ix = ilon0(ii);
%             fprintf(fid,'%14.5f%14.5f%12.4f%12.4f%12.4f\n',lon(ix)+360,lat(ix),ve(ix),vn(ix),vmag(ix) );   
%         end
%         fclose(fid);

        %-------------------------------------
        % compute euler poles (lon-lat) and anti-poles (lon-lat)
        [elat,elon] = xyz2latlon(exyz);
        [elat_anti,elon_anti] = antipode(elat,elon);

        % write euler poles to file (lon-lat)
        ww = [smod '_fix_' stref '_euler_lonlat.dat'];
        ofile = [dir_surface ww];
        disp(['writing to file ' ofile]);
        fid = fopen(ofile,'w');
        for ii=1:nump, fprintf(fid,'%14.5f %14.5f \n',elon(ii),elat(ii) ); end
        fclose(fid);

        % write euler anti-poles to file (lon-lat)
        ww = [smod '_fix_' stref '_euler_anti_lonlat.dat'];
        ofile = [dir_surface ww];
        disp(['writing to file ' ofile]);
        fid = fopen(ofile,'w');
        for ii=1:nump, fprintf(fid,'%14.5f %14.5f \n',elon_anti(ii),elat_anti(ii) ); end
        fclose(fid);
        %-------------------------------------

        % write the global velocity field for the plate ipick for fixed ifix
        if ipick_figs==1
            % KEY COMMAND: pick plate for which you compute the global
            % velocities (mm/yr) for the reference frames with ifix fixed
            %ipick = 11;  % 11 is PAC
            
            ipick_vec = [11 12 13 37 12];  % PAC
            %ipick_vec = [10 11 11 32 11];  % NAM
            ipick = ipick_vec(imodel);
            disp('  '); disp([' example plate is ' names{ipick}]);

            % compute global ipick surface velocity
            evec = exyz(:,ipick);
            Pxyz = latlon2xyz(lat,lon,earthr);
            Vrtp = euler2gps(evec, Pxyz);
            ve_uniform = Vrtp(3,:)';
            vn_uniform = -Vrtp(2,:)';
            vmag_uniform = sqrt( ve_uniform.^2 + vn_uniform.^2 );

            % for GMT plotting, set zero values to <0
            veps = 1e-4; vmag_uniform(vmag_uniform < veps) = -veps;

            if q <= 4
                % write plate model vector COMPONENTS (fixed ifix plate) to file (mm/yr)
                ww = [smod '_' stref '_' name_labs{ipick} 'only_vec_q' stq '.dat'];
                ofile = [dir_surface ww];
                disp(['writing to file ' ofile]);
                fid = fopen(ofile,'w');
                for ii = 2:num-1  % exclude the poles
                    fprintf(fid,'%14.5f %14.5f %12.4f %12.4f \n',...
                        lon(ii),lat(ii),ve_uniform(ii),vn_uniform(ii));   
                end
                fclose(fid);
            end

            % write plate model vector MAGNITUDES (fixed ifix plate) to file (mm/yr)
            ww = [smod '_' stref '_' name_labs{ipick} 'only_mag_q' stq '.dat'];
            ofile = [dir_surface ww];
            disp(['writing to file ' ofile]);
            fid = fopen(ofile,'w');
            for ii = 1:num
                fprintf(fid,'%14.5f%14.5f%12.4f \n',...
                    lon(ii),lat(ii),vmag_uniform(ii));   
            end
%             % repeat the lon=0 points at lon=360 for GMT plotting with grdimage
%             ilon0 = find( abs(lon) < 0.01);
%             for ii = 1:length(ilon0)
%                 ix = ilon0(ii);
%                 fprintf(fid,'%14.5f%14.5f%12.4f \n',...
%                     lon(ix)+360,lat(ix),vmag_uniform(ix) );   
%             end
%             fclose(fid);

        end
    end  % iwrite

end  % irow

%==========================================================================

if 0==1
    % By default, the program will figure out which plate each input point is on,
    % and then it will return the surface velocity. This next example shows
    % how one can get the surface velocity of a point that is outside the
    % boundary of the plate.
    clear, clc, close all
    % input points
    lon = [-144 -150];
    lat = [59 61];
    
    % choose:
    % a) plate model (with reference frame)
    % b) whether to fix a plate (this eliminates the relevance of the reference frame)
    % c) the plate for which you calculate surface velocities
    imodel = 10;
    [exyz,names,name_labs,dir_bounds,ssfx,smod] = get_plate_model(imodel,true);
    
    % four examples
    for iex = 1:4
        iplate = [];
        switch iex
            case 1, ifix = 13;              % fixed North America
            case 2, ifix = 99;              % default reference frame
            case 3, ifix = 99; iplate = 1;  % Pacific motion
            case 4, ifix = 99; iplate = 13; % North America motion
        end
        opts = {0,1,0,iplate};

        % compute velocity field for VECTORS (coarse mesh)
        % note: plots a vector field using the quiver command
        [lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] ...
            = platemodel2gps(lon,lat,imodel,ifix,opts);
        disp(sprintf('iex = %i',iex));
        for ii=1:length(lon)
            if isempty(iplate), ip = iplate_vec(ii); else ip = iplate; end
            disp(sprintf('input point (%7.2f, %7.2f) is (%.2f, %.2f) [%.2f mm/yr] on %s (%s) [%s, ifix=%i]',...
                lon(ii),lat(ii),ve(ii),vn(ii),sqrt(ve(ii)^2 + vn(ii)^2),names{ip},name_labs{ip},smod,ifix));
        end
        figure; quiver(lon,lat,ve,vn);
    end
    
% input point (-144.00,   59.00) is (-18.07, 52.25) [55.29 mm/yr] on Pacific (pa) [morvel_nnr, ifix=13]
% input point (-150.00,   61.00) is (0.00, -0.00) [0.00 mm/yr] on North_America (na) [morvel_nnr, ifix=13]
% 
% input point (-144.00,   59.00) is (-27.97, 31.58) [42.19 mm/yr] on Pacific (pa) [morvel_nnr, ifix=99]
% input point (-150.00,   61.00) is (-8.08, -21.64) [23.10 mm/yr] on North_America (na) [morvel_nnr, ifix=99]
%   
% input point (-144.00,   59.00) is (-27.97, 31.58) [42.19 mm/yr] on Pacific (pa) [morvel_nnr, ifix=99]
% input point (-150.00,   61.00) is (-28.82, 32.07) [43.12 mm/yr] on Pacific (pa) [morvel_nnr, ifix=99]
%   
% input point (-144.00,   59.00) is (-9.90, -20.67) [22.92 mm/yr] on North_America (na) [morvel_nnr, ifix=99]
% input point (-150.00,   61.00) is (-8.08, -21.64) [23.10 mm/yr] on North_America (na) [morvel_nnr, ifix=99]

    
end

%==========================================================================
