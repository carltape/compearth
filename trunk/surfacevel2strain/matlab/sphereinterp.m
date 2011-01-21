%
% moho_model.m
% Carl Tape, 17-Jan-2011
%
% Run script for sphereinterp_grid.m and sphereinterp_est.m.
%
% calls sphereinterp_grid.m, sphereinterp_est.m, utm2ll.m
% called by xxx
%

clc
clear
close all
format short, format compact

% TEMP
path(path,'/home/carltape/matlab/scripts');

earthr = 6371*1e3;      % earth radius (m)
deg = 180/pi;
msize = 6^2;

iwrite = 0;
iwavelet = 1;   % =1 for estimation; =0 to view data only
idata = 1;      % index for dataset

% 1-1 california moho
% 2-1 socal moho
% 1-2 USGS crystaline basement
% 1-4 SJB base Tertiary
% 3-3 Maricopa basement
% 4-5 nenana gravity
ropt  = input(' Type an index corresponding to a region (1=cal, 2=socal, 3=maricopa, 4=nenana): ');
dopt  = input(' Type an index corresponding to a dataset (1=moho, 2,3,4, 5=grav): ');

%========================================================
% GET 1D DATA SET

[dlon,dlat,d,dsig,ax0,slabel] = get_1D_dataset_carl(ropt,dopt);

%====================================================================
% ESTIMATE A SMOOTH MOHO MAP USING SPHERICAL WAVELETS

if iwavelet==1
    
    switch dopt
        case 1            
            qmin = 2; qmax = 8; % qmax = 8 or 9
            nlam = 40; ilampick = 2;
            ntrsh = 3;
            nx = 200;
        case 2           
            qmin = 2; qmax = 7;
            nlam = 40; ilampick = -10;  % hand-pick lambda
            ntrsh = 3;
            nx = 200;
            file0 = '/home/carltape/GOCAD/surfaces/cal_basement_outline.dat';
            [polyx,polyy,polylon,polylat] = textread(file0,'%f%f%f%f');
            
        case 3
            qmin = 5; qmax = 11;   % qmax = 11 or 12
            nlam = 40; ilampick = 1;
            ntrsh = 3;
            nx = 100;
            file0 = '/home/carltape/MOHO/DATA/SJB_gocad/SouthernBasement_poly.dat';
            [polyx,polyy] = textread(file0,'%f%f');
            [polylon,polylat] = utm2ll(polyx,polyy,szone,1);    
            
        case 4            
            qmin = 2; qmax = 10;   % qmax = 9 or 10
            nlam = 40; ilampick = 1;
            ntrsh = 3;
            nx = 200;
            file0 = '/home/carltape/MOHO/DATA/SJB_gocad/bt_poly.dat';
            [polyx,polyy] = textread(file0,'%f%f');
            [polylon,polylat] = utm2ll(polyx,polyy,szone,1);
            
        case 5
            % can take up to 30 minutes
            qmin = 5; qmax = 11;   % may want q=12 near the basin
            nlam = 40; ilampick = 1;
            ntrsh = 5;
            nx = 200;
    end
    
    nx = 50; qmin = 2; qmax = 7;   % testing
    qsec = round(mean([qmin qmax]));
    qparm = {qmin,qsec,qmax,ntrsh};
    rparm = {nlam,ilampick};
    if exist('polylon','var'), pparm = {nx,polylon,polylat}; else pparm = {nx}; end
    
    % KEY COMMAND: call sphereinterp.m
    %[dest,dest_plot,lam0,dlon_plot,dlat_plot,na,nb] = ...
    %    sphereinterp(dlon,dlat,d,dsig,ax0,qparm,rparm,pparm);
    
    [spline_tot] = sphereinterp_grid(dlon,dlat,ax0,qparm);
    
    [dest,dest_plot,lam0,dlon_plot,dlat_plot,na,nb] = ...
        sphereinterp_est(spline_tot,dlon,dlat,d,dsig,ax0,rparm,pparm);

end

%----------------------------------------------------------------
% WRITE FILES

if and(iwavelet==1,iwrite==1)
    %idata = 1;   % index for input data set
    %igrid = 1;   % index for ax0
    ftag = sprintf('%s_q%2.2i_q%2.2i_id%2.2i_ig%2.2i',slabel,qmin,qmax,idata,igrid);
    ndata = length(dest);
    nplot = length(dest_plot);
    
    dir_output = '/home/carltape/MOHO/WAVELET/MATLAB_EST/';
    %flab = [dir_output slabel '_' stqtag{1} '_' sprintf('ic%2.2i_im%2.2i',idata,sub_opt) ];
    flab = [dir_output ftag];
    
    disp('writing files with tag:'); disp(flab);
    
    % data and estimated field
    fid = fopen([flab '.dat'],'w');
    stfmt = '%18.8e%18.8e%18.8e%18.8e%18.8e\n';
    for ii=1:ndata
        fprintf(fid,stfmt,dlon(ii),dlat(ii),d(ii),dest(ii),dsig(ii));
    end
    fclose(fid);

    % estimated field for a regular grid
    fid = fopen([flab '_plot.dat'],'w');
    stfmt = '%18.8e%18.8e%18.8e\n';
    for ii=1:nplot
        fprintf(fid,stfmt,dlon_plot(ii),dlat_plot(ii),dest_plot(ii));
    end
    fclose(fid);

    % write bounds to file
    fid = fopen([flab '_bounds.dat'],'w');
    fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',ax0(1),ax0(2),ax0(3),ax0(4));
    fclose(fid);
    
    % write regularization parameter to file
    fid = fopen([flab '_lambda.dat'],'w');
    fprintf(fid,'%18.8e\n',lam0);
    fclose(fid);
end

%========================================================
