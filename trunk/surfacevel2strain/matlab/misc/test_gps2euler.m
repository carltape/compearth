%
% test_gps2euler.m
% Carl Tape, 21-August-2007
%
% This function tests gps2euler.m, which seeks to determine a best-fitting
% euler vector from a local velocity field on the sphere.
%
% An older version of this program is commented below.
%
% calls global2local.m
% called by xxx
%

clc
clear
close all
format short
format compact

% add path to additional matlab scripts (specify bdir)
user_path;

idata = 0;      % =1 to run the NASA REASON dataset
                % =0 to run on a synthetic dataset

deg = 180/pi;
earthr = 6371*1e3;    % earth radius, meters
                              
%-----------------------------------------
    
% KEY: input points
% The size of the region, the number of points, and the perturbation of the
% synthetic field from the actual field will all have a major effect on the
% ability to compute the best-fitting euler vector.
if idata==0
    n = 100;
    lon0 = randomvec(-122,-117,n);
    lat0 = randomvec(30,34,n);
    %lon0 = randomvec(-180,180,n);
    %lat0 = randomvec(-90,90,n);
    lon = lon0;
    lat = lat0;
else
    % NASA REASON GPS velocity field
    [name,lon0,lat0,start_date,ve0,vn0,vu0,se0,sn0,su0] = read_gps_REASON;
    sigma_horz = sqrt( se0.^2 + sn0.^2 );
    vmag_horz  = sqrt( ve0.^2 + vn0.^2 );
    percent_error_horz = 100 * sigma_horz ./ vmag_horz;
    
    % extract observations that will be used to find a best-fitting pole
    max_percent_error = 100;      % maximum percent error (horizontal only)
    max_west_lon = -109;        % maximum west position
    min_south_lon = 30;        % minimum south position
    ikeep = find( and( and( vn0 <= 0, percent_error_horz <= max_percent_error),...
        and( lon0 >= max_west_lon, lat0 >= min_south_lon)));
    
    % toss out observations in Colorado
    %inds = getsubset(lon0,lat0,[-106 -104 39 41]);
    %ikeep = setdiff(ikeep,inds);
    
    %disp([ikeep lon0(ikeep) lat0(ikeep)]);
    disp(' Datapoints for which we determine the best-fitting euler vector:');
    disp('     lon         lat       ve       vn       vmag      sigma    sigma-percent ');
    disp(sortrows([lon0(ikeep) lat0(ikeep) ve0(ikeep) vn0(ikeep) vmag_horz(ikeep) sigma_horz(ikeep) percent_error_horz(ikeep)],2));
    
    lon = lon0(ikeep);
    lat = lat0(ikeep);
    ve  = ve0(ikeep);
    vn  = vn0(ikeep);
    se  = se0(ikeep);
    sn  = sn0(ikeep);
    n   = length(ikeep);
    
    %se  = zeros(n,1);
    %sn  = zeros(n,1);
end

%Pxyz0 = latlon2xyz(lat0,lon0,earthr);
%Pxyz = latlon2xyz(lat,lon,earthr);

if idata==0
    Pxyz = latlon2xyz(lat,lon,earthr);
    
    % get a set of euler vectors for a plate model
    elatlon_in = [ -50.3837  107.8820    0.7554]'  % euler vector
    exyz_in = euler_convert(elatlon_in,0);
    %[lon, lat, ve, vn, iplate_vec, exyz, names, name_labs] = platemodel2gps(lon,lat,3,11,{0,0,0});
    %exyz = [-0.1479    0.4584   -0.5819]';
    %elatlon = euler_convert(exyz,1)
    Vrtp = euler2gps(exyz_in, Pxyz);      % surface vel (local r,th,ph), MM/YR
    vn0  = -Vrtp(2,:)';
    ve0  = Vrtp(3,:)';
    vmag0 = sqrt(ve0.^2 + vn0.^2);

    % add noise to the field
    npert = 0.1;
    nmag = randomvec(0,npert,length(lon)) .* vmag0;
    naz  = randomvec(0,360,length(lon));
    vecs = [ve0 vn0] + [nmag.*cos(naz/deg) nmag.*sin(naz/deg)];
    ve = vecs(:,1);
    vn = vecs(:,2);
    %disp([ve0 ve vn0 vn]);
    
    se = zeros(n,1);
    sn = zeros(n,1);
end

vmag = sqrt(ve.^2 + vn.^2);

figure; hold on;
if idata==0, quiver(lon,lat,ve0,vn0,'r'); end
quiver(lon,lat,ve,vn,'b');

figure; plot(vmag,'.');
xlabel(' observation index');
ylabel(' magnitude of horizontal velocity, mm/yr');

%----------------
% comput best-fitting euler vector, and then the predicted field

se = zeros(n,1);
sn = zeros(n,1);
[evec, ve_est, vn_est, Pxyz] = gps2euler(lon,lat,earthr,ve,vn,se,sn);
%elatlon_out = gps2euler(Pxyz,ve,vn,se,sn)       % KEY COMMAND
%exyz_out = euler_convert(elatlon_out,0)

break

% predictions at the points used in the inverse problem
Vrtp_out = euler2gps(exyz_out, Pxyz);
vn_out  = -Vrtp_out(2,:)';
ve_out  = Vrtp_out(3,:)';
ve_res = ve_out - ve;
vn_res = vn_out - vn;

figure; hold on;
quiver(lon,lat,ve,vn,'b');
quiver(lon,lat,ve_out,vn_out,'r');

figure; hold on;
quiver(lon,lat,ve_res,vn_res,'k');
title(' residuals');

%--------------------

if idata == 1
    % predictions at ALL points
    Vrtp0_out = euler2gps(exyz_out, Pxyz0);
    vn0_out  = -Vrtp0_out(2,:)';
    ve0_out  = Vrtp0_out(3,:)';
    ve0_res = ve0_out - ve0;
    vn0_res = vn0_out - vn0;

    %figure; hold on;
    %quiver(lon0,lat0,ve0,vn0,'b');
    %quiver(lon0,lat0,ve0_out,vn0_out,'r');

    figure; hold on;
    quiver(lon0,lat0,ve0_res,vn0_res,'k');
    title(' residuals');
end

%===================================================================
% OLD VERSION

%--------------------------------------------------
% PROBLEM: Given the local surface velocities v1,v2 at points r1,r2,
% determine w, the rotation vector used to compute v1,v2.
%
% SOLUTION #1: If v1 and v2 are linearly independent, then the direction of
% w is v1 x v2 (sign still undetermined), since w is orthogonal to both v1 and v2.
% Then you know sin t1, where t1 is the angle between w and r1, since you know
% the direction of w. Then from |v1|=|w x r1|=|w||r1|sin t1 you get |w|.
% Then you have w except for sign. Choose the sign so that determinant(w,r,v) > 0,
% but make sure that the coordinate system you are using when you get the
% coordinates of w, r, and v is itself right-handed. 
%
% SOLUTION #2: Use SOLUTION #1 to get an initial guess for the Euler pole,
% and then construct a least squares problem to find the best-fitting pole.
% THE PROBLEM IS LINEAR, BUT I OPTED TO USED THE GENERAL NON-LINEAR LEAST
% SQUARES PROGRAM genfit.m.
%--------------------------------------------------

% stfile = '/home/carltape/splines/';
% 
% %----------------------------------
% 
% ww2 = 'socal_gps_field'; load([stfile ww2 '.dat']);
% socal_field = eval(ww2);
% lon = socal_field(:,1);
% lat = socal_field(:,2);
% %vph = socal_field(:,3);
% %vth = -socal_field(:,4);
% 
% % pick an euler pole (fixed NAM), |E| = deg/Myr
% ePAC = [-0.2282 0.4057 -0.6882]';
% evec = ePAC * 1e-6 * pi/180;        % convert deg/Myr --> rad/yr
% [elat,elon] = xyz2latlon(evec);
% 
% ipick = 1;
% if ipick==1
%     % pick two lat-lon points
%     lon = lon(1:2);
%     lat = lat(1:2);
%     Pxyz = latlon2xyz(lat,lon,r);  % m
%     P1 = Pxyz(:,1);
%     P2 = Pxyz(:,2);
%     
%     % surface velocities, m/yr
%     V1 = cross(evec,Pxyz(:,1));
%     V2 = cross(evec,Pxyz(:,2));
%     
% elseif ipick==2
%     % find the points on the PAC plate (fix NAM)
%     [ve,vn,iPAC,iNAM] = gpsref_socal(lat, lon, [9 0 1 1]);
%     lat = lat(iPAC);
%     lon = lon(iPAC);
%     Pxyz = latlon2xyz(lat,lon,r);
%     
%     % find surface velocities, mm/yr (euler pole in deg/Myr) --
%     % pick the two largest vectors from the set
%     Vrtp = euler2gps(ePAC, Pxyz);
%     Vxyz = global2local(Vrtp,Pxyz,0);
%     [uV,vlen] = unit(Vxyz);
%     temp = sortrows([vlen' [1:length(vlen)]' ]);
%     %i1 = temp(end,2); i2 = temp(end-1,2);
%     i1 = temp(randomint(1,length(vlen),1),2); i2 = temp(randomint(1,length(vlen),1),2);
%     
%     V1 = Vxyz(:,i1) * 1e-3;     % m/yr
%     V2 = Vxyz(:,i2) * 1e-3;     % m/yr
%     P1 = Pxyz(:,i1);            % m
%     P2 = Pxyz(:,i2);            % m
% end
% 
% P1,P2,V1,V2
% norm(P1), norm(P2), norm(V1), norm(V2)
% 
% % surface velocities, m/yr
% V1rtp = global2local(V1, P1, 1);
% V2rtp = global2local(V2, P2, 1);
% 
% %-------------------------
% 
% % this gets the direction of E, which is in the direction V1 x V2
% [uE,Elen] = unit( cross(V1,V2) )
% 
% % this gets sin(t1), where t1 is the angle between uE and P1
% % |P1 x E| = |P1||uE|sin(t1) = |P1|sin(t1)
% norm(cross(P1,uE))
% sint1 = norm(cross(P1,uE)) / norm(P1)
% 
% % this gets |E| via |V1|=|E x P1|=|E||P1|sin t1 
% Emag = norm(V1) / (norm(P1)*sint1);
% E0 = uE * Emag;
% 
% % determine the sign of E0
% if det([ E0 P1 V1]) < 0, E0 = -E0; end
% 
% % convert: rad/yr --> deg/Myr
% Erec = E0 * 1e6 * 180/pi;
% 
% %-------------------------
% 
% stfm = '%.4f';
% st1 = num2str(sprintf(stfm, Erec(1)));
% st2 = num2str(sprintf(stfm, Erec(2)));
% st3 = num2str(sprintf(stfm, Erec(3)));
% st4 = num2str(sprintf(stfm, norm(Erec) ));
% stE = [' E ( x, y, z )  =  (' st1 ', ' st2 ', ' st3 '),  |E|  =  ' st4 ];
% 
% % using the computed euler pole, compute the surface velocities (m/yr)
% if 1==1
%     ZV1 = cross(E0,P1);
%     ZV2 = cross(E0,P2);
%     ZV1rtp = global2local(ZV1,P1,1);
%     ZV2rtp = global2local(ZV2,P2,1);
% else
%     ZVrtp = euler2gps(Erec, [P1 P2]) * 1e-3;
%     ZV1rtp = ZVrtp(:,1);
%     ZV2rtp = ZVrtp(:,2);
% end
% 
% % check with the original values:
% V1rtp - ZV1rtp
% V2rtp - ZV2rtp
% ePAC - Erec
% 
% % check
% [xlat,xlon] = xyz2latlon(Erec);
% [plat,plon] = xyz2latlon(P1);
% adist = arcdist(xlat,xlon,plat,plon);
% d = sin( adist/deg );
% d - sint1
% 
% % plot directions on unit sphere
% X = unit([P1 P2 V1 V2 Erec]);
% [Xlat, Xlon] = xyz2latlon(X);
% figure; nr=2; nc=1;
% subplot(nr,nc,1);
% globefun3(1,elat,elon,1,'b');
% globefun3(1,Xlat(1),Xlon(1),1,'k');
% globefun3(1,Xlat(2),Xlon(2),1,'k');
% globefun3(1,Xlat(5),Xlon(5),1,'r--');
% title({' P1, P2, E0, E',[' |E P1| = ' num2str(sprintf('%.2f', adist)) '^\circ' ], stE});
% subplot(nr,nc,2);
% globefun3(1,elat,elon,1,'b');
% globefun3(1,Xlat(3),Xlon(3),1,'k');
% globefun3(1,Xlat(4),Xlon(4),1,'k');
% globefun3(1,Xlat(5),Xlon(5),1,'r--');
% title(' V1, V2, E0, E');
% fontsize(10), orient tall, wysiwyg
% 
% disp('--------------------------------------');
% 
% % full velocity field
% lon0 = socal_field(:,1);
% lat0 = socal_field(:,2);
% ve0 = socal_field(:,3);
% vn0 = socal_field(:,4);
% Pxyz0 = latlon2xyz(lat0,lon0,r);
% num0 = length(lat0);
% 
% % get the LOCATIONS of the stations to get the Pacific Plate stations
% [ve,vn,iPAC,iNAM] = gpsref_socal(lat0,lon0,[9 0 1 1]);
% socal_PAC = socal_field(iPAC,:);
% 
% % get the largest N vectors on the Pacific Plate
% nbig = 40;
% vmag = sqrt( socal_PAC(:,3).^2 + socal_PAC(:,4).^2 );
% ibig = sortrows([vmag [1:length(vmag)]'], -1);
% igood = ibig(1:nbig,2);
% num = length(igood);
% stbig = [num2str(nbig) ' longest SoCal vectors'];
% 
% socal_PAC_big = socal_PAC(igood,:);
% lon  = socal_PAC_big(:,1);
% lat  = socal_PAC_big(:,2);
% ve   = socal_PAC_big(:,3);
% vn   = socal_PAC_big(:,4);
% Pxyz = latlon2xyz(lat,lon,r);
% Vrtp = [zeros(num,1) -vn ve]';
% 
% %-----------------------------------------
% % non-linear least squares fitting
% % THIS IS A LINEAR PROBLEM, BUT THIS WORKS FOR NOW.
% 
% % DATA
% data = [vn ; ve];
% 
% % INITIAL MODEL
% % initial guess euler pole = plate model euler pole
% m0 = ePAC;
% parm = Pxyz;
% V0 = feval('theory_euler2gps', m0, parm);
% vn_est0 = V0(1:num);
% ve_est0 = V0(num+1:end);
% 
% figure; hold on; s=0.01;
% quiver(lon,lat,s*ve,s*vn,0,'b');
% quiver(lon,lat,s*ve_est0,s*vn_est0,0,'r'); axis equal, grid on;
% xlabel(' Latitude'); ylabel(' Longitude'); 
% fontsize(10), orient tall, wysiwyg
% title({[' blue = data (' stbig ')'],' ',...
%     'red = initial estimated v-field based on global plate model (i.e, using ePAC)'});
% 
% % perturbations
% jogvec = 1e-6 *ones(1,3);
% 
% m1 = m0; itmx = 6;
% for ii=1:itmx
%     % KEY FUNCTION CALL:
%     [m1 e1] = genfit('theory_euler2gps', m1, jogvec, data, parm);
%     m1
% 
%     % error measurements
%     % rms residual = sqrt[ sum of squares of res / number of res ]
%     est = feval('theory_euler2gps', m1, parm);
%     res = data - est;
%     rms = sqrt( (res' * res) / length(res) );
%     stRMS = [' RMS = ' num2str(sprintf('%.4e', rms)) ';'];
%     disp(stRMS);
% end
% disp('Best-fit model:');
% m1
% disp(stRMS);
% 
% %-----------------------------------------
% 
% % estimated velocity at the nbig stations
% vn_est = est(1:num);
% ve_est = est(num+1:end);
% 
% % full velocity fields: data and model
% est_all = feval('theory_euler2gps', m1, Pxyz0);
% vn_est_all = est_all(1:num0);
% ve_est_all = est_all(num0+1:end);
% 
% figure; hold on; s=0.01;
% quiver(lon,lat,s*ve,s*vn,0,'b');
% quiver(lon,lat,s*ve_est,s*vn_est,0,'r'); axis equal, grid on;
% xlabel(' Latitude'); ylabel(' Longitude'); 
% title({[' blue = data (' stbig ')'],' ',...
%       'red = final estimated v-field after LSQ (i.e., using ePAC-fit)'});
% fontsize(10), orient tall, wysiwyg
% 
% % plot the full velocity fields: data, model, residual
% figure; hold on;
% quiver(lon0,lat0,s*ve0,s*vn0,0,'b');
% quiver(lon0,lat0,s*ve_est_all,s*vn_est_all,0,'r'); axis equal, grid on;
% xlabel(' Latitude'); ylabel(' Longitude');
% title([' blue = data,  red = model (fit to ' stbig ')']);
% fontsize(10), orient tall, wysiwyg
% 
% figure; hold on;
% quiver(lon0,lat0,ve0-ve_est_all,vn0-vn_est_all,'b');
% plot(lon0,lat0,'k.','markersize',10); axis equal;
% xlabel(' Latitude'); ylabel(' Longitude');
% title([' Residual  (' num2str(num) ' SoCal Pacific points should be approximately zero)']);
% fontsize(10), orient tall, wysiwyg
% 
% % plot the inital-guess and best-fit euler pole
% [X,Xlen] = unit([ePAC m1]);
% [Xlat, Xlon] = xyz2latlon(X);
% figure;
% globefun3(Xlen(1),Xlat(1),Xlon(1),1,'r');
% globefun3(Xlen(2),Xlat(2),Xlon(2),1,'b');
% title({[' red = ePAC (plate model Pacific euler pole),  |ePAC|  =  ' num2str(sprintf('%.3f', norm(ePAC))) ' deg/Myr'],' ',...
%        [' blue = ePAC-fit (SoCal-based Pacific euler pole),  |ePAC-fit|  =  ' num2str(sprintf('%.3f', norm(m1))) ' deg/Myr']});
% fontsize(10), orient tall, wysiwyg

%========================================================
