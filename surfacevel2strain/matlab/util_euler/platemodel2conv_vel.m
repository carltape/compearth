function [vvec_c, vmat_c, vmat_up, vmat_sub, gamma, exyz, names, name_labs] = ...
    platemodel2conv_vel(iups,isubs,az,lon,lat,imodel,ifix)
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
%   az      azimuth relative to which convergence angle is measured
%   lon     longitude of target point
%   lat     latitude of target point
%   imodel  index of plate model
%   ifix    index of fixed plate
%
% OUTPUT
%   vvec_c     angle relative to trench-normal: az_sub, az_up, az_conv = az_sub - az_up
%   vmat_c     n x 6: ve, vn, ve_normal, vn_normal, ve_parallel, vn_parallel
%   vmat_up    n x 6: ve, vn, ve_normal, vn_normal, ve_parallel, vn_parallel
%   vmat_sub   n x 6: ve, vn, ve_normal, vn_normal, ve_parallel, vn_parallel
%   gamma      obliquity angle, in degrees
%   exyz       plate model: euler vectors
%   names      plate model: plate names
%   name_labs  plate mode: plate labels
%
% Carl Tape, 16-May-2006
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
vvec_c   = NaN*ones(num,3);

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
        vaz_local = [0 -cos(az(ii)/deg) sin(az(ii)/deg)]';
        vaz_global = global2local(vaz_local, Pxyz, 0);
        
        % trench parallel vector (to the LEFT of the azimuth vector when
        % looking from the top)
        vpl_local = [0 -cos((az(ii)-90)/deg) sin((az(ii)-90)/deg)]';
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
        vvec_c(ii,1) = dot(vaz_global,Vxyz_sub);
        vvec_c(ii,2) = dot(vaz_global,Vxyz_up);
        vvec_c(ii,3) = vvec_c(ii,1) - vvec_c(ii,2);
        
        % obliquity from -180 to 180, measured from the azimuth vector
        % direction of azimuth vector is 0
        % coordinate system is E and N
        v1 = [sin(az(ii)/deg) cos(az(ii)/deg)];
        v2 = [vmat_sub(ii,1)-vmat_up(ii,1) vmat_sub(ii,2)-vmat_up(ii,2)];
        cosgam(ii) = dot(v1,v2)/(norm(v1)*norm(v2));
        csign(ii) = -sign(dot(vmat_c(ii,1:2), [sin((az(ii)-90)/deg) cos((az(ii)-90)/deg)]));
        
        %disp([ii v1 v2 dot(v1,v2) cosgam(ii) acos(cosgam(ii))*deg acos(cosgam(ii))*deg*csign(ii)]);
    end
end

% correct to get the right convention on obliquity: range is -180 to 180
gamma = csign .* acos(cosgam)*deg;

ifig = 1;
if ifig==1
    figure; nr=2; nc=1;
    subplot(nr,nc,1); scatter(wrapTo360(lon),lat,12^2,az,'filled'); caxis([0 360]); colorbar
    title('local trench-normal angle, degrees');
    subplot(nr,nc,2); scatter(wrapTo360(lon),lat,12^2,gamma,'filled'); caxis(60*[-1 1]); colorbar
    title('convergence angle, degrees');
end

%==========================================================================
% EXAMPLE

if 0==1
    % subset of Aleutian/Alaskan subduction transects from Lallemand2005
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
    az    = dall(:,3);
    lon   = dall(:,4);
    lat   = dall(:,5);
    
    imodel = 7;
    ifix = 99;
    [vvec_c, vmat_c, vmat_up, vmat_sub, gamma, exyz, names, name_labs] = ...
        platemodel2conv_vel(iups,isubs,az,lon,lat,imodel,ifix);
    for kk=1:length(az)
       disp(sprintf(' trench-normal = %6.1f, sub = %6.1f, up %6.1f, conv = %6.1f, oblique = %6.1f',...
           az(kk),vvec_c(kk,:),gamma(kk) ));
    end
end

%==========================================================================
