function [vvec_c, vmat_c, vmat_up, vmat_sub, gamma, exyz, names, name_labs] = ...
    platemodel2conv_vel(iups,isubs,azref,lon,lat,imodel,ifix)
%PLATEMODEL2CONV_VEL compute convergent velocities between plates
%
% This function computes the convergence velocity at a plate boundary,
% given the two plates at the boundary and the azimuth perpendicular to the
% boundary (pointing in the direction of the subducting plate in the case
% of a subduction boundary).
%
% INPUT
%   iups    index of upper plate
%   isubs   index of lower plate
%   azref   azimuth relative to which the convergence angle is measured
%   lon     longitude of target point
%   lat     latitude of target point
%   imodel  index of plate model
%   ifix    index of fixed plate
%
% OUTPUT
%   vvec_c     n x 6: angles (SEE CONVENTIONS BELOW)
%   vmat_c     n x 6: ve, vn, ve_normal, vn_normal, ve_parallel, vn_parallel
%   vmat_up    n x 6: ve, vn, ve_normal, vn_normal, ve_parallel, vn_parallel
%   vmat_sub   n x 6: ve, vn, ve_normal, vn_normal, ve_parallel, vn_parallel
%   gamma      obliquity angle, in degrees
%   exyz       plate model: euler vectors
%   names      plate model: plate names
%   name_labs  plate mode: plate labels
%
% Carl Tape, 2006-05-16
%

% constants
deg = 180/pi;
earthr = 6371*1e3;

% might need to make sure that longitude is [-180,180]
lat = lat(:);
lon = lon(:);
num = length(lat);

% get euler vectors for plate model
get_plate_model;

% display info on euler poles
if 1==1
    
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

%========================================================
% COMPUTE SURFACE VELOCITY FIELD

disp('  ');
if ifix==99
    disp(' leave v-field in the original reference frame');
    flab = ' in original reference frame';
else
    % fix one of the plates by subtracting its euler vector from all vectors
    disp([' fixed plate is ' name_labs{ifix} ' (' names{ifix} ')']);
    flab = [' with fixed ' name_labs{ifix}];
    exyz0 = exyz;
    exyz = exyz0 - repmat( exyz0(:, ifix), 1, nump );
end

% initialization
% 6 columns: ve, vn, ve_normal, vn_normal, ve_parallel, vn_parallel
vmat_up  = NaN*ones(num,6);     % upper plate
vmat_sub = NaN*ones(num,6);     % subducting plate
vmat_c   = NaN*ones(num,6);     % convergent velocity

% trench normal convergent velocity
vvec_c   = NaN*ones(num,6);

% obliquity
cosgam = NaN*ones(num,1);
gamma = NaN*ones(num,1);
csign = NaN*ones(num,1);

for ii=1:num
    % index of upper plate and subducting plate
    iup = iups(ii);
    isub = isubs(ii);
    
    if ~isnan([iup isub])
        % euler vectors for subducting plate and upper plate
        esub = exyz(:,isub);
        eup  = exyz(:,iup);
        
        % input point
        Pxyz = latlon2xyz(lat(ii),lon(ii),earthr);

        % azimuth: trench normal vector (in the direction of the subducting plate)
        vaz_local = [0 -cos(azref(ii)/deg) sin(azref(ii)/deg)]';
        vaz_global = global2local(vaz_local, Pxyz, 0);
        
        % trench parallel vector
        % (to the LEFT of the azimuth vector when looking from the top)
        vpl_local = [0 -cos((azref(ii)-90)/deg) sin((azref(ii)-90)/deg)]';
        vpl_global = global2local(vpl_local, Pxyz, 0);  
        
        %--------------
        % subducting plate
        
        % compute surface velocity field -- Vsub
        [Vrtp, Vxyz_sub] = euler2gps(esub, Pxyz);       % surface vel (local r,th,ph)
        vmat_sub(ii,2) = -Vrtp(2,:)';
        vmat_sub(ii,1) = Vrtp(3,:)';
        
        % determine projection vectors: trench normal
        temp = global2local(dot(vaz_global,Vxyz_sub)*unit(vaz_global),Pxyz,1);
        vmat_sub(ii,3:4) = [temp(3) -temp(2)];
        
        % determine projection vectors: trench parallel
        temp = global2local(dot(vpl_global,Vxyz_sub)*unit(vpl_global),Pxyz,1);
        vmat_sub(ii,5:6) = [temp(3) -temp(2)];
        
        %--------------
        % upper plate
        
        % compute surface velocity field -- Vup
        [Vrtp, Vxyz_up] = euler2gps(eup, Pxyz);         % surface vel (local r,th,ph)
        vmat_up(ii,2) = -Vrtp(2,:)';
        vmat_up(ii,1) = Vrtp(3,:)';
        
        % determine projection vectors: trench normal
        temp = global2local(dot(vaz_global,Vxyz_up)*unit(vaz_global),Pxyz,1);
        vmat_up(ii,3:4) = [temp(3) -temp(2)];
        
        % determine projection vectors: trench parallel
        temp = global2local(dot(vpl_global,Vxyz_up)*unit(vpl_global),Pxyz,1);
        vmat_up(ii,5:6) = [temp(3) -temp(2)];
        
        %--------------
        % convergent velocity (Vc = Vsub - Vup)
        
        vmat_c(ii,1:2) = vmat_sub(ii,1:2) - vmat_up(ii,1:2);    % vector
        vmat_c(ii,3:4) = vmat_sub(ii,3:4) - vmat_up(ii,3:4);    % trench normal
        vmat_c(ii,5:6) = vmat_sub(ii,5:6) - vmat_up(ii,5:6);    % trench parallel
        
        % trench-normal convergence velocity (note sign convention)
        vvec_c(ii,1) = dot(vaz_global,Vxyz_sub);        % ALWAYS POSITIVE
        vvec_c(ii,2) = dot(vaz_global,Vxyz_up);         % ALWAYS POSITIVE
        vvec_c(ii,3) = vvec_c(ii,1) - vvec_c(ii,2);
        
        % obliquity from -180 to 180, measured from the azimuth vector
        % direction of azimuth vector is 0
        % coordinate system is E and N
        v1 = [sin(azref(ii)/deg) cos(azref(ii)/deg)];
        v2 = [vmat_sub(ii,1)-vmat_up(ii,1) vmat_sub(ii,2)-vmat_up(ii,2)];
        cosgam(ii) = dot(v1,v2)/(norm(v1)*norm(v2));
        csign(ii) = -sign(dot(vmat_c(ii,1:2), [sin((azref(ii)-90)/deg) cos((azref(ii)-90)/deg)]));
        
        %disp([ii v1 v2 dot(v1,v2) cosgam(ii) acos(cosgam(ii))*deg acos(cosgam(ii))*deg*csign(ii)]);
    end
end

% correct to get the right convention on obliquity: range is -180 to 180
gamma = csign .* acos(cosgam)*deg;

% calculate azimuths of all vectors
ph_c = cart2pol(vmat_c(:,1),vmat_c(:,2));
ph_c = ph_c*180/pi;
az_c = ph2az(ph_c);
ph_up = cart2pol(vmat_up(:,1),vmat_up(:,2));
ph_up = ph_up*180/pi;
az_up = ph2az(ph_up);
ph_sub = cart2pol(vmat_sub(:,1),vmat_sub(:,2));
ph_sub = ph_sub*180/pi;
az_sub = ph2az(ph_sub);

vvec_c(:,4) = az_c;
vvec_c(:,5) = az_up;
vvec_c(:,6) = az_sub;

bfigure_basic = true;
bfigure_all = false;

if bfigure_basic
    if num > 1
        if max(lon) - min(lon) > 180, plon = wrapTo360(lon); else plon = lon; end
        figure; hold on;
        quiver(plon,lat,vmat_c(:,1),vmat_c(:,2),0,'k');
        quiver(plon,lat,vmat_up(:,1),vmat_up(:,2),0,'r');
        quiver(plon,lat,vmat_sub(:,1),vmat_sub(:,2),0,'b');

        figure; nr=2; nc=1;
        subplot(nr,nc,1); scatter(plon,lat,12^2,azref,'filled'); caxis([0 360]); colorbar
        title('local trench-normal angle, degrees');
        subplot(nr,nc,2); scatter(plon,lat,12^2,gamma,'filled'); caxis(60*[-1 1]); colorbar
        title('convergence angle, degrees');
    end
    
    % all on one figure
    figure; hold on;
    for ii=1:num
       %annotation('arrow',[0 0],[vmat_c(ii,1) vmat_c(ii,2)],'Color','k');
       h1 = plot([0 vmat_c(ii,1)],[0 vmat_c(ii,2)],'k');
       h2 = plot([0 vmat_up(ii,1)],[0 vmat_up(ii,2)],'r');
       h3 = plot([0 vmat_sub(ii,1)],[0 vmat_sub(ii,2)],'b');
       plot([vmat_sub(ii,1) vmat_c(ii,1)],[vmat_sub(ii,2) vmat_c(ii,2)],'r--');
       axis equal; grid on;
       xlabel('East velocity, mm/yr');
       ylabel('North velocity, mm/yr');
       legend([h1 h2 h3],'convergence','upper','subduct');
       if num==1
           title(sprintf('lon %.2f lat %.2f azcon %.1f azsub %.1f azsub %.1f [azref %.1f] %.1f %.1f %.1f',...
             lon(ii),lat(ii),vvec_c(ii,4:6),azref(ii),vvec_c(ii,1:3)));
       end
       if num > 40, continue; end
    end
end

if bfigure_all
    % individual figures
    for ii=1:num
       figure; hold on;
       plot([0 vmat_c(ii,1)],[0 vmat_c(ii,2)],'k');
       plot([0 vmat_up(ii,1)],[0 vmat_up(ii,2)],'r');
       plot([0 vmat_sub(ii,1)],[0 vmat_sub(ii,2)],'b');
       plot([vmat_sub(ii,1) vmat_c(ii,1)],[vmat_sub(ii,2) vmat_c(ii,2)],'r--');
       axis equal; grid on;
       xlabel('East velocity, mm/yr');
       ylabel('North velocity, mm/yr');
       legend([h1 h2 h3],'convergence','upper','subduct');
       title(sprintf('lon %.2f lat %.2f azcon %.1f azsub %.1f azsub %.1f [azref %.1f] %.1f %.1f %.1f',...
           lon(ii),lat(ii),vvec_c(ii,4:6),azref(ii),vvec_c(ii,1:3)));
       if num > 40, continue; end
    end
end

%==========================================================================
% EXAMPLE

if 0==1
    %% subset of Aleutian/Alaskan subduction transects from Lallemand2005
    clear, clc, close all
    dall = [
       32.0000   37.0000   28.0000  173.0000   51.6000
       32.0000   37.0000   16.0000  175.0000   51.1000
       32.0000   37.0000   21.0000  177.0000   50.8000
       32.0000   37.0000   12.0000  179.0000   50.4000
       32.0000   37.0000  351.0000 -179.0000   50.3000
       32.0000   37.0000  350.0000 -177.0000   50.4000
       32.0000   37.0000  352.0000 -175.0000   50.5000
       32.0000   37.0000  353.0000 -173.0000   50.9000
       32.0000   37.0000  342.0000 -171.0000   51.1000
       32.0000   37.0000  336.0000 -169.0000   51.5000
       32.0000   37.0000  338.0000 -167.0000   52.0000
       32.0000   37.0000  338.0000 -165.0000   52.5000
       32.0000   37.0000  340.0000 -163.0000   53.1000
       32.0000   37.0000  344.0000 -161.0000   53.5000
       32.0000   37.0000  342.0000 -159.0000   53.8000
       32.0000   37.0000  341.0000 -157.0000   54.2000
       32.0000   37.0000  323.0000 -155.0000   54.8000
       32.0000   37.0000  330.0000 -153.0000   55.6000
       32.0000   37.0000  331.0000 -151.0000   56.2000
       32.0000   37.0000  320.0000 -149.0000   57.1000
       32.0000   37.0000  315.0000 -147.0000   58.0000
       32.0000   37.0000  311.0000 -145.5000   59.1000
       32.0000   37.0000  345.0000 -144.5000   59.2000
   ];
    iups  = dall(:,1);
    isubs = dall(:,2);
    azref = dall(:,3);
    lon   = dall(:,4);
    lat   = dall(:,5);
    
    %azref = azref*0;  % TEST THE CASE OF IGNORING THE TRENCH ORIENTATION
    imodel = 7;
    ifix = 99;
    [vvec_c, vmat_c, vmat_up, vmat_sub, gamma, exyz, names, name_labs] = ...
        platemodel2conv_vel(iups,isubs,azref,lon,lat,imodel,ifix);
    for kk=1:length(azref)
       disp(sprintf(' trench-normal = %6.1f, sub = %6.1f, up %6.1f, conv = %6.1f, oblique = %6.1f',...
           azref(kk),vvec_c(kk,:),gamma(kk) ));
    end
    
    %% Cook Inlet and Susitna region
    iups = 32;      % NA
    isubs = 37;     % PA
    azref = 0;
    lon = -151.5;
    lat = 61.5;
    imodel = 7;
    ifix = 99;
    [vvec_c, vmat_c, vmat_up, vmat_sub, gamma, exyz, names, name_labs] = ...
        platemodel2conv_vel(iups,isubs,azref,lon,lat,imodel,ifix);
    
    
end

%==========================================================================
