%
% function get_gps_dataset.m
% CARL TAPE, 25-March-2009
% printed xxx
%
% This function loads a velocity field for certain points on the sphere.
% 
% calls platemodel2gps.m, read_gps_3D.m
% called by surfacevel2strain.m
%

function [dlon,dlat,vu,vs,ve,su,sn,se,ax1,dir,slabel,stref] ...
    = get_gps_dataset(ropt,dopt,istore,iplate_model)

% directories
dir0       = '/home/carltape/SURFACEVEL2STRAIN/';
dir_gps    = [dir0 'gps_data/'];
dir_grids  = [dir0 'fortran/grids_output/'];

% geographic regions: lat-lon boxes specified in get_subgrids.f90
% THESE GRIDPOINTS SHOULD BE GENERATED IN get_subgrids.f90
% irow is only relevant when dealing with a plate model
% --> add dopt index to the full list of possibilities below
ropts_all = {'west_us','cal','socal','taiwan','tibet','cascadia','asia','parkfield','japan','wedge'};
irow_all = [11 11 11 17 7 7 11 11 11 11];
slabel = ropts_all{ropt};
irow = irow_all(ropt);

nropt = length(ropts_all);
if length( find(ropt == [1:nropt]) )==0, error(' check region options (ropt)'); end
if length( find(dopt == [0 1 2 3 4 10:13 20:23 30:33 40:43 50:53 60:63 70:73 80:83]) )==0
    error(' check data options (dopt)');
end
if dopt == 0, istore = 0; end
sdopt = sprintf('d%2.2i', dopt);

% % file label and region (see test_platemodel2gps.m)
% % THESE GRIDPOINTS SHOULD BE GENERATED IN get_subgrids.f90
% % irow is only relevant when dealing with a plate model
% switch ropt
%     case 1, slabel = 'west_us';  irow = 11;
%     case 2, slabel = 'cal';      irow = 11;
%     case 3, slabel = 'socal';    irow = 11;
%     case 4, slabel = 'taiwan';   irow = 17;
%     case 5, slabel = 'tibet';    irow = 7;
%     case 6, slabel = 'cascadia'; irow = 7;
% end

% directory containing the bounds and the spline gridpoints
sub_opt = 1;   % index in get_subgrids.f90
dir = [dir_grids 'subgrids_' slabel '_' sprintf('%2.2i', sub_opt) '/'];

% load the bounds for the region
ax1 = load([dir 'bounds_input.dat']);
lonmin = ax1(1); lonmax = ax1(2);
latmin = ax1(3); latmax = ax1(4);

if istore == 1   % use specific v-field data (velocities in MM/YR)

    % OBSERVED VELOCITY FIELDS
    if dopt == 1        % NASA REASON dataset (modified in REASON_gps_dat.m)
        %filename = [dir_gps 'US/reason_fixed_NAM_subset_3D.dat'];
        %filename = [dir_gps 'US/reason_fixed_NAM_v2_subset_3D.dat'];
        filename = [dir_gps 'US/reason_subset_3D.dat'];

    elseif dopt == 2    % CCMM, v1.0 (modified in socal_gps_dat.m)
        filename = [dir_gps 'US/california/socal_vfield_4p0_3D.dat'];
        
    elseif dopt == 3    % Jean-Phillipe, central Asia
        filename = [dir_gps 'ASIA/asia/data_JP.txt'];
        
        [name,dlon,dlat,vn,ve,sn,se] = textread(filename,'%s%f%f%f%f%f%f');
        dlon = lonshift(dlon,[1 1]);
        vs = -vn;
        ndata = length(dlon);
        vu = zeros(ndata,1);
        su = zeros(ndata,1); sn = zeros(ndata,1); se = zeros(ndata,1);
        
    elseif dopt == 4    % Takeo Ito, Japan
        filename = [dir_gps 'ASIA/japan/japan_takeo_ito_subset_3D.dat'];

    % SYNTHETIC VELOCITY FIELDS
    elseif dopt >= 10
        filename = [dir_gps 'synthetic/syn_vfield_' sdopt '_3D.dat'];
    end
    
    % NOTE: It is simpler to store the velocity field in a format that can
    % be read by read_gps_3D.m
    disp([' about to read ' filename]);
    [dlon,dlat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] = read_gps_3D(filename);
    
    ndata = length(dlon);

    % combine into one array
    data_all = [dlon dlat vu vn ve su sn se];

    % ELIMINATE OBSERVATIONS OUTSIDE THE SPECIFIED REGION
    i_inregion = getsubset(dlon,dlat,ax1);
    data_all = data_all(i_inregion,:);

    % convert to meters (per second)
    data_all(:,[3:8]) = data_all(:,[3:8])*1e-3;

    dlon = data_all(:,1);
    dlat = data_all(:,2);
    vu   = data_all(:,3); vs = -data_all(:,4); ve = data_all(:,5);
    su   = data_all(:,6); sn =  data_all(:,7); se = data_all(:,8);

    if 0==1
        % see socal_gps_syn.m to add Gaussian noise to the velocity field
    end
    
    % dummy variable to send back
    stref = 'NA';

else   % not using GPS data, so get plate model velocity field

    % generate uniform mesh for PLATE MODEL fields
    numx = 60;
    [dlon,dlat] = gridvec(lonmin,lonmax,numx,latmin,latmax);

    % you may want to constrain datapoints OUTSIDE the 'target region'
    %dfac = 0.5 * (latmax-latmin);
    %[dlon2,dlat2] = gridvec(lonmin-dfac,lonmax+dfac,numx,latmin-dfac,latmax+dfac);
    
    %irow = 11;  % FIXED PLATE

    ifix_mat = [ 1 1 1 2        %  1 AFR
                 2 2 4 6        %  2 ANT
                 3 3 5 7        %  3 ARA
                 4 4 6 8        %  4 AUS
                 5 5 7 13       %  5 CAR
                 6 6 99 15      %  6 COC
                 7 7 8 18       %  7 EUR
                 99 8 9 21      %  8 IND
                 8 9 99 22      %  9 JDF
                 9 10 11 32     % 10 NAZ
                 10 11 11 32    % 11 NAM
                 11 12 13 37    % 12 PAC
                 12 13 14 39    % 13 PHI (PS)
                 13 15 16 46    % 14 SAM
                 99 14 99 42    % 15 SCO
                 99 99 17 48    % (16) SU
                 99 99 18 52    % (17) YA
                 ];
    ifix = ifix_mat(irow, iplate_model);
    opts = {0,1,0};

    % compute velocity field for VECTORS
    [dlon, dlat, ve, vn, iplate_vec, exyz, names, name_labs] ...
        = platemodel2gps(dlon,dlat,iplate_model,ifix,opts);

    % convert to from mm/yr to m/yr
    ve = ve * 1e-3;               % ve = v_phi
    vs = -vn * 1e-3;              % vs = v_theta

    stref = name_labs{ifix};      % fixed plate

    % initialize other vectors
    ndata = length(dlon);
    vu = zeros(ndata,1);
    su = zeros(ndata,1); sn = zeros(ndata,1); se = zeros(ndata,1);
end

%===================================================================
