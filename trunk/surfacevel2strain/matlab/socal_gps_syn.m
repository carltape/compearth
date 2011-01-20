%
% socal_gps_syn.m
% Carl Tape, 12-June-2009
%
% This program generates a synthetic GPS velocity field based on a simple
% analytical model of strike-slip faulting.
%
% Several of these fields were analyzed in 
%   Tape, Muse, Simons, Dong, Webb, "Multiscale Estimation of GPS velocity
%   fields," Geophysicsl Journal International, 2009.
%
% SYNTHETIC VELOCITY FIELDS (isyn_field):
%   10 -- strike-slip, uniform field, no errors
%   11 -- strike-slip, uniform field, errors
%   12 -- strike-slip, irregular field, no errors
%   13 -- strike-slip, irregular field, errors
%   20-23 -- rotational
%   30-33 -- 3D finite strike-slip
%   40-43 -- 3D infinite strike-slip
%   50-53 -- 3D coseismic subduction event
%   60-63 -- microplate rotation, I
%   70-73 -- microplate rotation, II
%   80-83 -- volcanic source (Mogi)
%
% NOTE: the UTM conversion function requires the Matlab Mapping Toolbox.
%
% calls
%    strike_slip_sphere.m
%    platemodel2gps.m
%    get_gps_dataset.m
%    socal_gps_split.m
% called by xxx
%

clear
close all

% add path to additional matlab scripts
path(path,[pwd '/util']);
path(path,[pwd '/func']);
path(path,[pwd '/misc/okada']);

deg = 180/pi;
earthr = 6371*1e3;
izone = 11;
szone = '11S';

bdir = '/home/carltape/compearth/surfacevel2strain/';   % USER change

iwrite = 0;

%===========================================
% USER INPUT
ropt = input(' Enter an index corresponding to a region (1=us, 2=cal, 3=socal, ..., 9=japan, etc): ');
dopt = input(' Enter an index that controls the type of synthetic field (1-80): ');

dmat = NaN*ones(60,3);
dmat(10,:) = [1 2 0]; dmat(11,:) = [1 2 1]; dmat(12,:) = [1 1 0]; dmat(13,:) = [1 1 1];
dmat(20,:) = [2 2 0]; dmat(21,:) = [2 2 1]; dmat(22,:) = [2 1 0]; dmat(23,:) = [2 1 1];
dmat(30,:) = [3 2 0]; dmat(31,:) = [3 2 1]; dmat(32,:) = [3 1 0]; dmat(33,:) = [3 1 1];
dmat(40,:) = [4 2 0]; dmat(41,:) = [4 2 1]; dmat(42,:) = [4 1 0]; dmat(43,:) = [4 1 1];
dmat(50,:) = [5 2 0]; dmat(51,:) = [5 2 1]; dmat(52,:) = [5 1 0]; dmat(53,:) = [5 1 1];
dmat(60,:) = [6 2 0]; dmat(61,:) = [6 2 1]; dmat(62,:) = [6 1 0]; dmat(63,:) = [6 1 1];
dmat(70,:) = [7 2 0]; dmat(71,:) = [7 2 1]; dmat(72,:) = [7 1 0]; dmat(73,:) = [7 1 1];
dmat(80,:) = [8 2 0]; dmat(81,:) = [8 2 1]; dmat(82,:) = [8 1 0]; dmat(83,:) = [8 1 1];

isyn_field  = dmat(dopt,1);  % see list above
iobs_points = dmat(dopt,2);  % 1=NASA REASON; 2=uniform
ierrors     = dmat(dopt,3);  % 1=add errors; 0=no errors

%iobs_points = input(' Enter a set of observations points (1 = NASA REASON; 2 = uniform): ');
%isyn_field = input(' Enter a type of synthetic velocity field (1 = strike-slip fault; 2 = rotational field): ');
%ierrors = input(' Enter 1 to include errors (0 otherwise) : ');

if any(isnan([isyn_field iobs_points ierrors])==1), error(' NaN value for booleans'); end

% number of points for plotting synthetic fault (great circle)
ngc = 1000;

%-------------------------------------------

% load the SAF for plotting: isaf, latsaf, lonsaf, xsaf, ysaf
[lonsaf,latsaf,xsay,ysaf] = textread('../gmt/input/safdata2.dat','%f%f%f%f');
%load('/home/carltape/matlab/scripts/safdata2');
nsaf = length(lonsaf);

%===========================================

% To get a realistic set of observation points, we load a real GPS dataset.
% NOTE: units in METERS and METERS/YR
dopt_0 = 1;   % NASA REASON cGPS dataset (408 pts in socal)
%dopt_0 = 2;   % CCMMv1 (1093 pts in socal)
[lon_gps,lat_gps,vu_gps,vs_gps,ve_gps,su,sn,se,ax0,slabel,stref] ...
    = get_gps_dataset(ropt,dopt_0,1,0);

lonmin = ax0(1); lonmax = ax0(2);
latmin = ax0(3); latmax = ax0(4);
ax1 = ax0;

% labels for files
sdopt = sprintf('%2.2i',dopt);
file1 = ['syn_vfield_d' sdopt];
%file1 = ['syn_vfield_' slabel '_d' sdopt];

% generate uniform mesh for SYNTHETIC velocity fields
if iobs_points == 2
    numx = 35;
    [lon_gps,lat_gps] = gridvec(lonmin,lonmax,numx,latmin,latmax);
end

n = length(lon_gps);
disp('  '); disp(file1);
disp([num2str(n) ' observation points']);

figure; plot(lon_gps, lat_gps, '.');

%---------------------------------------

if isyn_field == 1      % strike-slip fault

    dep = 15*1e3;       % locking depth, m
    slip = 0.035;        % slip rate, m/yr
    %lat1 = 37.25; lon1 = -122; lat2 = 32; lon2 = -114.75;
    lat1 = 36.0; lon1 = -120.7; lat2 = 35.6; lon2 = -120.0;
    [Vtheta,Vphi,xdat,Glat,Glon,evecA,evecB] = ...
        strike_slip_sphere(lat_gps,lon_gps,dep,slip,lat1,lon1,lat2,lon2);

    % rotation rate
    disp('rotation rates of the two hemispheres, rad/yr:');
    norm(evecA)*1e-6/deg
    norm(evecB)*1e-6/deg
    euler_convert(evecA,1)    %  41.5308, 9.3523, 0.1574
    euler_convert(evecB,1)    % -41.5308, -170.6477, 0.1574
    
    % no vertical component
    vu_an = zeros(n,1);
    vs_an = Vtheta;
    ve_an = Vphi;
    
    % index of the point to FIX for the zero-value of the velocity field
    % --> take point farthest from the SAF, and let this represent 'rigid North America'
    %[junk, imax] = min(xdat);

    inds = getsubset(Glon,Glat,ax0);
    lat_gc = Glat(inds);
    lon_gc = Glon(inds);
    
    stL = [' Strike-slip fault locked at z = ' num2str(dep/1e3) ' km'];
    stslip = [' Slip rate = ' num2str(slip*1e3) ' mm/yr'];
    stit = [stL ', ' stslip];
    
    figure; hold on;
    quiver(lon_gps,lat_gps,Vphi,-Vtheta); axis equal
    plot(Glon,Glat,'r');
    axis(ax0);
    
%     % pick starting and finishing point for the GREAT-CIRCLE fault,
%     % indexed from NW to SE
%     %istart = 1; iend = nsaf;
%     istart = 8; iend = 35;
% 
%     [latsaf_gc, lonsaf_gc] = gcwaypts(latsaf(istart),lonsaf(istart),latsaf(iend),lonsaf(iend));
%     az_start = azimuth(latsaf(istart),lonsaf(istart),latsaf(iend),lonsaf(iend));
% 
%     % pick the starting point of the fault to be some distance NW
%     dinc = 20;
%     [lattrk,lontrk] = track1(latsaf(istart),lonsaf(istart),az_start+180,dinc,'degrees');
%     lon_start = lontrk(iend);
%     lat_start = lattrk(iend);
%     
%     % compute new azimuth
%     az_start = azimuth(lat_start,lon_start,latsaf(iend),lonsaf(iend));
% 
%     %---------------------------------------
%     % fault parameters
% 
%     imodel = 1;
%     mtag = sprintf('%2.2i', imodel);
% 
%     % five-parameter model:
%     %  1, locking depth (km)
%     %  2, slip rate (cm/yr)
%     %  3-4-5, great circle of plane: (lon_start, lat_start, azimuth)
%     parm_mat = [15 -4.0 lon_start lat_start az_start];
%     parms = parm_mat(imodel,:);
% 
%     % SAF: 15 km, 3.5 cm/yr
% 
%     %---------------------------------------
% 
%     % top and bottom depth of fault (m)
%     Linc = Inf*1e3;     % fault length (infinite)
%     Ltop = parms(1)*1e3;      % depth of top of fault
%     Lbot = Ltop+Linc;   % depth of bottom of fault
% 
%     % analytical solution (see fault_2D.m, especially for axes convention)
%     dip         = 90;
%     depth       = Ltop;
%     W           = (Lbot-Ltop)/sin(dip/deg);
%     slip        = parms(2)*1e-2;        % + = right lateral, - = left lateral
%     fault_type  = 1;                % strike-slip
% 
%     %---------------------------------------
% 
%     % simply use a line of points
%     x0 = 2e5;
%     xdat = linspace(-x0,x0,1000);
% 
%     stx = 'x , horizontal distance from strike-slip fault (m)';
% 
%     [dh,dv] = fault_2D(xdat,slip,W,dip,depth,fault_type);
%     uyAN = dh;
%     uzAN = dv;
% 
%     % compute e_xy : note the axis conventions
%     %exyAN = gradient(uyAN*1e-2,xdat*1e3);       % all in meters
%     exyAN = 0.5*gradient(uyAN,xdat);
%     exyAN = -exyAN;
% 
%     %stL = {' Strike-slip fault extending from', [ num2str(Ltop/1e3) ' km to ' num2str(Lbot/1e3) ' km depth']};
%     stL = [' Strike-slip fault locked at z = ' num2str(Ltop/1e3) ' km'];
%     stslip = [' Slip rate v_y = ' num2str(slip*1e2) ' cm/yr'];
% 
%     figure; nr=3; nc=1;
%     xmin = min(xdat); xmax = max(xdat);
%     ymax = 1.2*max([ max(abs(uyAN)) max(abs(uzAN)) ]);
% 
%     stit = [stL ', ' stslip];
% 
%     clim = [-1 4]*1e-7;
% 
%     subplot(nr,nc,1); hold on;
%     plot(xdat,uyAN,'r.');
%     %plot([xmin xmax],[0 0],'k:');
%     %plot([0 0],[-ymax ymax],'k:');
%     axis([xmin xmax -ymax ymax]); grid on;
%     xlabel(stx); ylabel(' v_y ,  horizontal velocity (m/yr)');
%     set(gca,'YDir','reverse'); title(stit);
% 
%     subplot(nr,nc,2); hold on;
%     plot(xdat,uzAN,'r.');
%     axis([xmin xmax -ymax ymax]); grid on;
%     xlabel(stx); ylabel(' v_z ,  vertical velocity (m/yr)'); title(stit);
% 
%     subplot(nr,nc,3); hold on;
%     plot(xdat,exyAN,'r.');
%     axis([xmin xmax -1e-7 4e-7]); grid on;
%     xlabel(stx); ylabel(' Strain rate, \epsilon_{xy}'); title(stit);
% 
%     fontsize(9), orient tall; wysiwyg
% 
%     %=========================
% 
%     %------------
%     % get the distance-to-fault for all the points
%     % get the points corresponding to GPS station locations
% 
%     figure; hold on; axis equal;
%     plot(lon_gps, lat_gps, 'k.');
% 
%     % GC defining the strike of the PLANAR San Andreas Fault (see above)
%     lon_start = parms(3); lat_start = parms(4); az_start = parms(5);
%     rng = 25;
%     [lat_gc,lon_gc] = track1(lat_start,lon_start,az_start,rng,[],'degrees',ngc);
% 
%     % AZIMUTH varies along the great circle
%     az_gc = zeros(ngc,1);
%     for ii=1:ngc-1
%         az_gc(ii) = azimuth(lat_gc(ii),lon_gc(ii),lat_gc(ii+1),lon_gc(ii+1));
%     end
%     az_gc(ngc) = az_gc(ngc-1);
% 
%     plot(lon_gc, lat_gc, 'g--');
%     plot(lonsaf, latsaf, 'r.-', 'linewidth', 2);
%     plot(lonsaf_gc, latsaf_gc, 'r--');
% 
%     % apply finite rotation (1) to the great circle; (2) to the points
%     elatlon = [lat_start lon_start az_start+180];
%     [latgc_rot, longc_rot, R] = euler_rot_tec(lat_gc,lon_gc,elatlon);
%     [latgps_rot, longps_rot, R] = euler_rot_tec(lat_gps,lon_gps,elatlon);
% 
%     plot(longps_rot, latgps_rot, 'k+');
%     plot(longc_rot, latgc_rot, 'r--');
% 
%     % compute the distance between (1) a point on the rotated great circle that
%     % WAS at the same latitude as the original station and (2) the now-rotated station
%     gcrot_lon = mean(longc_rot)*ones(n,1);
%     gcrot_lat = latgps_rot;
%     [dist_deg,az_deg] = distance(gcrot_lat, gcrot_lon, latgps_rot, longps_rot);
% 
%     % find the index of the point on the rotated GC that is closest to each
%     % observation point
%     i_az = zeros(n,1);
%     for ii=1:n
%         [dtemp,i_az(ii)] = min(distance(latgps_rot(ii),longps_rot(ii),latgc_rot,longc_rot));
%     end
%     az_gps = az_gc(i_az);    % from NW to SE
% 
%     % convert distances to km, and get the sign, based on the azimuth
%     dist = earthr/deg * dist_deg;
%     iSW = (az_deg > 180);
%     dist(iSW) = -dist(iSW);
% 
%     %dist = arcdist(gcrot_lat, gcrot_lon, latgps_rot, longps_rot) ...
%     %    / deg*earthr .* sign(longps_rot - gcrot_lon);
% 
%     figure; hold on;
%     [X,Y,Z] = griddataXB(lon_gps,lat_gps,dist,100,'cubic');
%     pcolor(X,Y,Z); shading interp;
%     plot(lon_gc,lat_gc,'k','linewidth',2);
%     colorbar;
% 
%     xdat = dist;
%     %------------
% 
%     [dh,dv] = fault_2D(xdat,slip,W,dip,depth,fault_type);
%     uyAN = dh;
%     uzAN = dv;
% 
%     figure; nr=2; nc=1;
%     xmin = min(xdat); xmax = max(xdat);
%     ymax = 1.2*max([ max(abs(uyAN)) max(abs(uzAN)) ]);
% 
%     stit = [stL ', ' stslip];
% 
%     clim = [-1 4]*1e-7;
% 
%     subplot(nr,nc,1); hold on;
%     plot(xdat,uyAN,'r.');
%     %plot([xmin xmax],[0 0],'k:');
%     %plot([0 0],[-ymax ymax],'k:');
%     axis([xmin xmax -ymax ymax]); grid on;
%     xlabel(stx); ylabel(' v_y ,  horizontal velocity (m/yr)');
%     set(gca,'YDir','reverse'); title(stit);
% 
%     subplot(nr,nc,2); hold on;
%     plot(xdat,uzAN,'r.');
%     axis([xmin xmax -ymax ymax]); grid on;
%     xlabel(stx); ylabel(' v_z ,  vertical velocity (m/yr)'); title(stit);
% 
%     fontsize(9), orient tall; wysiwyg
% 
%     %----------------------
%     % KEY: to plot the vectors, we use the (+/-) magnitudes and the azimuth of the
%     % planar, vertical San Andreas Fault
% 
%     ve_an = uyAN .* cos((az_gps+90)/deg);
%     vs_an = uyAN .* sin((az_gps+90)/deg);
% 
%     % index of the point to FIX for the zero-value of the velocity field
%     % --> take point farthest from the SAF, and let this represent 'rigid North America'
%     [junk, imax] = max(xdat);
%     ve_ref = ve_an(imax);
%     vs_ref = vs_an(imax);
%     
%     % no vertical component
%     vu_an = zeros(n,1);

elseif isyn_field == 2      % rotational field

    iplate_model = 3;      % REVEL
    ifix = 11;             % fix NA plate
    [lon_gps, lat_gps, ve, vn, iplate_vec, exyz, names, name_labs] ...
        = platemodel2gps(lon_gps,lat_gps,iplate_model,ifix,{0,1,0});

    % we want a UNIFORM v-field, so just use the pacific (w.r.t NA) euler pole
    ipac = 13;
    evec1 = exyz(:,ipac);                  % euler vector (deg/Myr)
    evec2 = euler_convert(evec1,1);        % euler vector (wx,wy,wz)

    if 1==1
        % we MOVE the pole to within the SoCal region to create a highly
        % rotational velocity field; the -omega is to have the same sense
        % of rotation as the plate motion field
        omega = evec2(3);
        evec2_mod = [35 -116 -omega]';
        %evec2_mod = [35 -116 omega]';
        evec1_mod = euler_convert(evec2_mod,0);
        evec1 = evec1_mod;
        evec2 = evec2_mod;
    end
    stE = sprintf(' Euler vector: (lat = %.2f, lon = %.2f, omega = %.2f deg/Myr)',evec2(1),evec2(2),evec2(3));
    stit = stE;

    Pxyz = latlon2xyz(lat_gps,lon_gps,earthr);  % input points
    Vrtp = euler2gps(evec1, Pxyz);
    vs_an = Vrtp(2,:)' * 1e-3;                  % mm --> m
    ve_an = Vrtp(3,:)' * 1e-3;                  % mm --> m
    vmag  = sqrt(ve_an.^2 + vs_an.^2);

    % euler vector in radians per year
    evec_radpyr = evec1 * 1e-6 * pi/180;
    omega_radpyr = norm(evec_radpyr);       % 1.3184e-08

    figure; hold on;
    [X,Y,Z] = griddataXB(lon_gps, lat_gps, vmag, 100, 'cubic');
    pcolor(X,Y,Z); shading interp;
    quiver(lon_gps,lat_gps,ve_an,-vs_an,'k');
    plot(lonsaf, latsaf, 'k-','linewidth',1);
    axis equal; axis(ax0);
    colorbar
    xlabel(' Longitude'); ylabel(' Latitude');
    title(' Velocity field');
    %title(' Plate-motion velocity field (PA w.r.t fixed NA)');
    
    % no vertical component
    vu_an = zeros(n,1);

 elseif (isyn_field == 3 ||  isyn_field == 4)  % strike-slip fault

    % fault length, meters
    if isyn_field == 3
        L = 400*1e3;      % finite fault
    else
        L = 10000*1e3;    % infinite fault
    end

    % FAULT PARAMETERS (to generate 3D displacement field)
    strike = 120;
    dip = 90;
    rake = 0;
    sind = sin(dip/deg);
    W = 20*1e3;                     % width of fault plane, meters
    Ztop = 4*1e3;                   % depth of top of locked fault, meters
    Z = W*sind/2 + Ztop;            % depth of center of fault, meters
    Zbot = W*sind + Ztop;           % depth of bottom of locked fault, meters
    %Z = 30000;                     % depth of center of fault, meters
    
    slip3D = -3;           % fault slip, meters
    lambda = 33*1e9;        % Lame coefficient
    mu = 33*1e9;            % Lame coefficient
    slipelev = 0;
    nL = 1;
    nW = 1;
    
    stit = sprintf('W = %.0f km, L = %.0f km, Ztop = %.0f km, Zbot = %.0f km, S = %.1f m',...
        W*1e-3,L*1e-3,Ztop*1e-3,Zbot*1e-3,slip3D);
    
    % the okada codes asks for coordinates in utm (meters), in a local reference frame
    %izone = 11;
    szone = '11S';
    [gps_e,gps_n,i_zone] = utm2ll(lon_gps,lat_gps,szone,0);
    min_e = min(gps_e); max_e = max(gps_e);
    min_n = min(gps_n); max_n = max(gps_n);
    if 0==1
        center_e = min_e + 0.55*(max_e - min_e);
        center_n = min_n + 0.45*(max_n - min_n);
    else
        center_e = 4.428880869195670e+05;
        center_n = 3.783107758856259e+06;
    end
    
    gps_up = zeros(size(gps_e));
    tmpfile = 'misc/slip_constant.txt';
    rakeslipchiIguales(rake,slip3D,slipelev,nL,nW,tmpfile);
    [vn_an,ve_an,vu_an,ux,uy,uz,tmpvar] = okadaparchesNEV(lambda,mu,L,nL,W,nW,...
        Z, strike, dip, tmpfile, gps_e-center_e, gps_n-center_n, gps_up);
    
    % compute fault segment at surface
    utm_pt1 = [center_e center_n] + L/2*[cos((90-strike)/deg) sin((90-strike)/deg)];
    utm_pt2 = [center_e center_n] + L/2*[cos((90-strike+180)/deg) sin((90-strike+180)/deg)];
    [lon_pt1,lat_pt1] = utm2ll(utm_pt1(1,1),utm_pt1(1,2),szone,1);
    [lon_pt2,lat_pt2] = utm2ll(utm_pt2(1,1),utm_pt2(1,2),szone,1);
    [lat_gc,lon_gc] = gcwaypts(lat_pt1,lon_pt1,lat_pt2,lon_pt2,ngc);
    
    vs_an = -vn_an;   
    
elseif isyn_field == 5  % thrust fault

    % FAULT PARAMETERS (to generate 3D displacement field)
    strike = -5;
    dip = 15;
    rake = 90;
    sind = sin(dip/deg);
    L = 800*1e3;               % length of fault, meters
    W = 40*1e3;                % width of fault plane, meters
    Ztop = 5*1e3;              % depth of top of locked fault, meters
    Z = W*sind/2 + Ztop;       % depth of center of fault, meters
    Zbot = W*sind + Ztop;      % depth of bottom of locked fault, meters
    %Z = 30000;                % depth of center of fault, meters
    
    slip3D = 8;           % fault slip, meters
    lambda = 33*1e9;        % Lame coefficient
    mu = 33*1e9;            % Lame coefficient
    slipelev = 0;
    nL = 1;
    nW = 1;
    
    stit = sprintf('W = %.0f km, L = %.0f km, Ztop = %.0f km, Zbot = %.0f km, S = %.1f m',...
        W*1e-3,L*1e-3,Ztop*1e-3,Zbot*1e-3,slip3D);
    
    % the okada codes asks for coordinates in utm (meters), in a local reference frame
    %izone = 10;
    [gps_e,gps_n,i_zone] = utm2ll(lon_gps,lat_gps,szone,0);
    min_e = min(gps_e); max_e = max(gps_e);
    min_n = min(gps_n); max_n = max(gps_n);
    if 0==1
        center_e = min_e + 0.55*(max_e - min_e);
        center_n = min_n + 0.45*(max_n - min_n);
    else
        center_e = 3.681497582726048e+05 - 100*1e3;
        center_n = 4.907065202496171e+06 + 25*1e3;
    end
    
    gps_up = zeros(size(gps_e));
    tmpfile = 'slip_constant.txt';
    rakeslipchiIguales(rake,slip3D,slipelev,nL,nW,tmpfile);
    [vn_an,ve_an,vu_an,ux,uy,uz,tmpvar] = okadaparchesNEV(lambda,mu,L,nL,W,nW,...
        Z, strike, dip, tmpfile, gps_e-center_e, gps_n-center_n, gps_up);
    
    % compute fault segment at surface
    utm_pt1 = [center_e center_n] + L/2*[cos((90-strike)/deg) sin((90-strike)/deg)];
    utm_pt2 = [center_e center_n] + L/2*[cos((90-strike+180)/deg) sin((90-strike+180)/deg)];
    [lon_pt1,lat_pt1] = utm2ll(utm_pt1(1,1),utm_pt1(1,2),szone,1);
    [lon_pt2,lat_pt2] = utm2ll(utm_pt2(1,1),utm_pt2(1,2),szone,1);
    [lat_gc,lon_gc] = gcwaypts(lat_pt1,lon_pt1,lat_pt2,lon_pt2,ngc);
    
    vs_an = -vn_an;
    
elseif or(isyn_field == 6,isyn_field == 7)      % microplate rotation field
    
    vs_an = zeros(n,1);
    ve_an = zeros(n,1);
    
    % get points west and east of the San Andreas
    [iwest, ieast] = socal_gps_split(ax0,lon_gps,lat_gps);
    
    % specify euler vectors
    if isyn_field == 6
        evec_east = [36.0 -117.5 -0.5]';
        evec_west = [33.5 -119.0 1.0]';
    else
        evec_east = [36.0 -117.5 0.5]';
        evec_west = [33.5 -119.0 1.0]';
    end
    evecW = euler_convert(evec_west,2);
    evecE = euler_convert(evec_east,2);
    evec_radpyr_west = evecW * 1e-6 * pi/180;
    evec_radpyr_east = evecE * 1e-6 * pi/180;
    omega_radpyr_west = norm(evec_radpyr_west);
    omega_radpyr_east = norm(evec_radpyr_east);
    
    %stW = sprintf(' Euler vector: (lat = %.2f, lon = %.2f, omega = %.2f deg/yr)',evec_west(1),evec_west(2),evec_west(3));
    %stE = sprintf(' Euler vector: (lat = %.2f, lon = %.2f, omega = %.2f deg/yr)',evec_east(1),evec_east(2),evec_east(3));
    stit = sprintf(' Euler vectors: (%.2f, %.2f, %.2f deg/Myr) and (%.2f, %.2f, %.2f deg/Myr) ',...
        evec_west(1),evec_west(2),evec_west(3),evec_east(1),evec_east(2),evec_east(3));

    % Pacific field
    Pxyz = latlon2xyz(lat_gps(iwest),lon_gps(iwest),earthr);  % input points
    Vrtp = euler2gps(evecW, Pxyz);
    vs_an(iwest) = Vrtp(2,:)' * 1e-3;	% mm --> m
    ve_an(iwest) = Vrtp(3,:)' * 1e-3;	% mm --> m
    
    % North America field
    Pxyz = latlon2xyz(lat_gps(ieast),lon_gps(ieast),earthr);  % input points
    Vrtp = euler2gps(evecE, Pxyz);
    vs_an(ieast) = Vrtp(2,:)' * 1e-3;	% mm --> m
    ve_an(ieast) = Vrtp(3,:)' * 1e-3;	% mm --> m
    
    vmag  = sqrt(ve_an.^2 + vs_an.^2);

    figure; hold on;
    [X,Y,Z] = griddataXB(lon_gps, lat_gps, vmag, 100, 'cubic');
    pcolor(X,Y,Z); shading interp;
    quiver(lon_gps,lat_gps,ve_an,-vs_an,'k');
    plot(lonsaf, latsaf, 'k-','linewidth',1);
    axis equal; axis(ax0);
    colorbar
    xlabel(' Longitude'); ylabel(' Latitude');
    title(' Plate-motion velocity field (PA w.r.t fixed NA)');
    
    % no vertical component
    vu_an = zeros(n,1);
    
    % for GMT, plot the boundary between the microplates
    lon_gc = lonsaf;
    lat_gc = latsaf;
    
elseif isyn_field == 8         % volcanic inflation
    
    vs_an = zeros(n,1);
    ve_an = zeros(n,1);
    vu_an = zeros(n,1);
    
    %[xvec,yvec] = gridvec(-10,80,100,-10,75);
    %[Ux, Uy, Uz] = forward_mogi([10 15 -2 1], xvec, yvec);
    %figure; hold on; quiver(xvec,yvec,Ux,Uy,'k');
    %axis equal; title(' Volanic inflation field');
    
    % two volcanic signals
    lon0 = -116.5; lat0   = 34.5; depth0 = 100*1e3; c0 = 1e9;
    lon1 = -118.0; lat1   = 34.0; depth1 = 30*1e3; c1 = -1e7;
    %izone = 11;
    m0 = [lon0 lat0 depth0 c0]';
    m1 = [lon1 lat1 depth1 c1]';
    stvol0 = sprintf('(lon = %.1f, lat = %.1f, depth = %.1f km)',lon0,lat0,depth0*1e-3);
    stvol1 = sprintf('(lon = %.1f, lat = %.1f, depth = %.1f km)',lon1,lat1,depth1*1e-3);
    stit = ['Volcanic sources at ' stvol0 ' and ' stvol1];
    
    % compute displacement field at the input gridpoints
    %izone = 11;
    [ve_an0, vn_an0, vu_an0] = forward_mogi(m0, lon_gps, lat_gps, szone);
    [ve_an1, vn_an1, vu_an1] = forward_mogi(m1, lon_gps, lat_gps, szone);
    ve_an = ve_an0 + ve_an1;
    vu_an = vu_an0 + vu_an1;
    vs_an = -(vn_an0 + vn_an1);
    
    vmag  = sqrt(ve_an.^2 + vs_an.^2);

    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on;
    [X,Y,Z] = griddataXB(lon_gps, lat_gps, vmag, 100, 'cubic');
    pcolor(X,Y,Z); shading interp;
    quiver(lon_gps,lat_gps,ve_an,-vs_an,'k');
    axis equal; axis(ax0); colorbar
    xlabel(' Longitude'); ylabel(' Latitude');
    title({'Horizontal Field :',stit});

    subplot(nr,nc,2); hold on;
    [X,Y,Z] = griddataXB(lon_gps, lat_gps, vu_an, 100, 'cubic');
    pcolor(X,Y,Z); shading interp;
    axis equal; axis(ax0); colorbar
    xlabel(' Longitude'); ylabel(' Latitude');
    title({'Vertical Field :',stit});
    
end

if (isyn_field >= 3) && (isyn_field <= 5)
   plot_okada_fields(gps_e-center_e,gps_n-center_n,lambda,mu,L,nL,W,nW,Z,strike,dip,rake,slip3D,slipelev);
end

%======================================================================

figure; nr=2; nc=2;
subplot(nr,nc,1); plot(1e3*vu_an,'.'); title(' syn field: Vu (Vr), mm');
subplot(nr,nc,2); plot(1e3*vs_an,'.'); title(' syn field: Vs (Vtheta), mm');
subplot(nr,nc,3); plot(1e3*ve_an,'.'); title(' syn field: Ve (Vphi), mm');
subplot(nr,nc,4); plot(1e3*sqrt(vs_an.^2 + ve_an.^2),'.'); title(' syn field: Ve (Vphi), mm');

figure; quiver(lon_gps,lat_gps,ve_an,-vs_an); title(' synthetic field');

%----------------------
% ADD GAUSSIAN ERRORS IN METERS

if ierrors == 0
    % Even if you do not want errors, it is mathematically consistent with
    % least-squares to have Gaussian errors, and methods like
    % cross-validation are more effective with them (even if they're tiny).
    %sigma_val = 1e-6 * median(abs([vs_an ; ve_an]));
    sigma_val = 0;
    sigma_obs = sigma_val * ones(n,1);
    
else

    median(se)*1e3
    median(sn)*1e3
    median(su)*1e3
    
    % artificial erros in data, sigma in METERS
    if isyn_field==1                        % 2D interseismic strike-slip
        sigma_val = 0.5*1e-3;
        sigma_obs = sigma_val * ones(n,1);
        
    elseif isyn_field==2                    % rotational field
        %sigma_frac = 0.04;
        %sigma_obs = vmag * sigma_frac;
        sigma_val = 0.05*1e-3;
        sigma_obs = sigma_val * ones(n,1);
        
    elseif or(isyn_field==3, isyn_field==4) % 3D coseismic strike-slip
        sigma_val = 5*1e-3;
        sigma_obs = sigma_val * ones(n,1);
        
    elseif isyn_field==5                    % 3D coseismic thrust
        sigma_val = 5*1e-3;
        sigma_obs = sigma_val * ones(n,1); 
        
    elseif or(isyn_field==6, isyn_field==7)  % microplates rotation
        %sigma_frac = 0.04;
        %sigma_obs = vmag * sigma_frac; 
        sigma_val = 0.05*1e-3;
        sigma_obs = sigma_val * ones(n,1);
        
    elseif isyn_field==8                    % volcanic dilation
        sigma_val = 0.05*1e-3;
        sigma_obs = sigma_val * ones(n,1);
        
    end
end

sigma_mean = mean(sigma_obs);

% note randn --> normal distribution
vu_std = randn(n,1) .* sigma_obs;
vs_std = randn(n,1) .* sigma_obs;
ve_std = randn(n,1) .* sigma_obs;

if sigma_obs > 0
    figure; nr=3; nc=1;
    subplot(nr,nc,1); plot_histo(vu_std,[-4*sigma_mean: sigma_mean/2 : 4*sigma_mean]); ylabel(' Vu-error (m)');    
    subplot(nr,nc,2); plot_histo(vs_std,[-4*sigma_mean: sigma_mean/2 : 4*sigma_mean]); ylabel(' Vs-error (m)');
    subplot(nr,nc,3); plot_histo(ve_std,[-4*sigma_mean: sigma_mean/2 : 4*sigma_mean]); xlabel(' Ve-error (m)');
    orient tall; wysiwyg

    figure; quiver(lon_gps,lat_gps,ve_std,-vs_std); axis equal;
    title(' Errors to add to the synthetic velocity field');
    orient tall; wysiwyg
end

% add synthetic errors
vu_an = vu_an + vu_std;
vs_an = vs_an + vs_std;
ve_an = ve_an + ve_std;

figure; quiver(lon_gps,lat_gps,ve_an,-vs_an);
title(' synthetic field, referenced, with errors added');

%----------------------
% REMOVE A ROTATIONAL FIELD

% THIS IS NOT ADEQUATE ON THE SPHERE
%ve_an = ve_an - ve_ref;
%vs_an = vs_an - vs_ref;

figure; quiver(lon_gps,lat_gps,ve_an,-vs_an); title(' synthetic field, referenced');

%----------------------

%xmax = 5e5; xmin = -xmax;
%ymin = 3.3e6; ymax = ymin + xmax-xmin;
%ax1 = [xmin xmax ymin ymax];
%ax1 = [-122 -114 31 37];
%ax1 = [-126 -113.3 30 43.5];

figure; nr=2; nc=1;

if and(iobs_points == 1, isyn_field == 1)
    subplot(nr,nc,1); hold on;
    quiver(lon_gps,lat_gps,ve_gps,-vs_gps);
    plot(lonsaf, latsaf, 'r-','linewidth',2);
    axis equal; axis(ax1);
    title('GPS velocity data');
    xlabel(' West-East (m)'); ylabel(' South-North (m)');
end

subplot(nr,nc,2); hold on;
quiver(lon_gps,lat_gps,ve_an,-vs_an);
plot(lonsaf, latsaf, 'r-','linewidth',1);
if any(isyn_field == [1 3 5]), plot(lon_gc, lat_gc, 'r--','linewidth',2); end
%plot(lon_gps(imax), lat_gps(imax), 'ro');
%plot([0 0],[ymin ymax],'r--');
axis equal; axis(ax1);
title({'Synthetic velocity field (socal-gps-syn.m)',' ',stit});
xlabel(' West-East (m)'); ylabel(' South-North (m)');

% subplot(nr,nc,4); hold on;
% plot(xmat,ymat,'k+','markersize',3);
% plot(E,N,'r.');
% xlabel(' West-East (m)'); ylabel(' South-North (m)');
% title({'Mesh and stations: ',stpm});
% axis equal; axis(ax1);

fontsize(9), orient tall; wysiwyg

%-----------------------

if and(iobs_points == 1, isyn_field == 1)
    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot( ve_gps, ve_an, '.'); axis equal, grid on;
    subplot(nr,nc,2); plot( vs_gps, vs_an, '.'); axis equal, grid on;

    pscale = 1;  % uniform scaling of vectors for plotting

    figure; hold on;
    quiver(lon_gps,lat_gps,ve_gps,-vs_gps,pscale,'b');
    quiver(lon_gps,lat_gps,ve_an,-vs_an,pscale,'r');
    legend(' SCEC model','synthetic great-circle fault','location','southwest');
    plot(lonsaf, latsaf, 'k-','linewidth',1);
    plot(lon_gc, lat_gc, 'k--','linewidth',2);
    %plot(lon_gps(imax), lat_gps(imax), 'ro');
    %plot([0 0],[ymin ymax],'r--');
    axis equal; axis(ax1);
    title({'Crustal velocities for an infinite strike-slip fault (socal-gps-syn.m)',' ',stit});
    xlabel(' West-East (m)'); ylabel(' South-North (m)');
    fontsize(9), orient tall; wysiwyg
end

if iwrite == 0
    disp('Not writing any synthetic velocity fields.');
else
    odir = [bdir 'data/synthetic/'];   % USER CHANGE
    filetag = [odir file1];
    disp(['write files of the form ' filetag]);
    
    % convert to mm
    vu = 1e3*vu_an;
    ve = 1e3*ve_an;
    vn = -1e3*vs_an;
    sn = 1e3*sigma_obs;
    se = 1e3*sigma_obs;
    su = 1e3*sigma_obs;
    
    ren = zeros(n,1); reu = zeros(n,1); rnu = zeros(n,1);
    start_date = zeros(n,1); finish_date = zeros(n,1);
    name = []; for ii=1:n, name{ii} = sprintf('%4.4i',ii); end; name = name(:);

    write_xy_points(filetag,lon_gps,lat_gps);
    write_gps_psvelo(filetag,lon_gps,lat_gps,ve,vn,se,sn,ren,name);
    %write_gps_psxy_vert(filetag,fscale,lon_gps,lat_gps,vu,su);
    write_gps_3D(filetag,lon_gps,lat_gps,ve,vn,vu,...
        se,sn,su,ren,reu,rnu,start_date,finish_date,name);

    % great circle (segment) defining the synthetic fault
    if any(isyn_field == [1 3 5 6 7])
        ww = ['gps_gc_d' sdopt '.dat'];
        fid = fopen([odir ww],'w');
        for ii = 1:length(lon_gc)
            fprintf(fid,'%18.8f%18.8f\n',lon_gc(ii),lat_gc(ii));
        end
        fclose(fid);
    end
end

if 0==1
    [X,Y,Z] = griddataXB(lon_gps,lat_gps,vu,100,'linear');
    figure; hold on;
    pcolor(X,Y,Z); plot(lon_gps,lat_gps,'.'); colorbar
end

%============================================================
