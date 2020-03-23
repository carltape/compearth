%
% sphereinterp_carl.m
% 
% Carl's version of sphereinterp.m -- which has his example.
%
% calls get_1D_dataset.m, sphereinterp_grid.m, sphereinterp_est.m
% called by sphereinterp.m
%

clc, clear, close all
format short, format compact

% add path to additional matlab scripts (specify bdir)
user_path;
    
%========================================================
% USER PARAMETERS

iwavelet = 1;   % =1 for estimation; =0 to view data only
iwrite = 0;

% CARL EXAMPLES (ropt-dopt)
% 2-1 california moho
% 1-1 socal moho
% 1-2 USGS crystalline basement
% 1-4 SJB base Tertiary
% 3-3 Maricopa basement (socal)
% 4-5 nenana gravity
% 5-6 alaska moho
% 2-7 California rock data
% 8-8 Alaska shear-wave splitting
ropt  = input(' Type an index corresponding to a region (1=socal, 2=cal, 3=maricopa, 4=nenana, 5=alaska, 8=sak): ');
dopt  = input(' Type an index corresponding to a dataset (1=moho,2,3,4,5=grav,8=split): ');

%========================================================
% GET DATA SET

[dlon,dlat,d,dsig,ax0,slabel,ulabel,isource,ftags,d2,stout] = get_1D_dataset_carl(ropt,dopt);
dir_output = '/home/carltape/MOHO/WAVELET/MATLAB_EST/';

%====================================================================
% ESTIMATE A SMOOTH MOHO MAP USING SPHERICAL WAVELETS

if iwavelet==1
    
    % NOTE: Only option dopt=1 is available as the example
    switch dopt
        case 1            
            qmin = 2; qmax = 8; % qmax = 8 or 9
            nlam = 40; ilampick = 2;
            ntrsh = 3;
            nx = 100;
            
        case 2           
            qmin = 2; qmax = 7;
            nlam = 40; ilampick = -10;  % hand-pick lambda
            ntrsh = 3;
            nx = 200;
            file0 = '/home/carltape/GOCAD/surfaces/cal_basement_outline.dat';
            [polyx,polyy,polylon,polylat] = textread(file0,'%f%f%f%f');
            
        case 3
            szone = '11S';
            qmin = 5; qmax = 11;   % qmax = 11 or 12
            %qmax = 8;       % testing
            nlam = 40; ilampick = 1;
            ntrsh = 3;
            nx = 200;
            file0 = '/home/carltape/MOHO/DATA/SJB_gocad/SouthernBasement_poly.dat';
            [polyx,polyy] = textread(file0,'%f%f');
            [polylon,polylat] = utm2ll(polyx,polyy,szone,1);    
            
        case 4      
            szone = '11S';
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
            
        case 6
            qmin = 4; qmax = 8;   % may want q=12 near the basin
            nlam = 40; ilampick = 2;
            ntrsh = 3;
            nx = 200;
            
        case 7
            qmin = 0; qmax = 8;
            nlam = 40; ilampick = 2;
            ntrsh = 3;
            nx = 50;
            
        case 8
            % alaska splitting
            qmin = 6; qmax = 8;
            nlam = 40; ilampick = 2;
            ntrsh = 3;
            nx = 20;    % number of points controlling plotting grid
            
        case 10
            % note: the default choices will NOT adequately represent the
            % highest MMI regions -- to do this, you need to choose larger
            % qmax and lower value of lambda for regularization.
            % Or if you use low q values, you may want to decrease
            % minlampwr in sphereinterp_est.m
            qmin = 0; qmax = 7;
            nlam = 40;
            ilampick = 2;
            ntrsh = 3;
            nx = 200;
    end
    
    %nx = 50; qmin = 2; qmax = 7;   % testing
    
    qsec = round(mean([qmin qmax]));
    qparm = {qmin,qsec,qmax,ntrsh};
    rparm = {nlam,ilampick};
    if exist('polylon','var')
        pparm = {nx,ulabel,polylon,polylat};
    else
        pparm = {nx,ulabel};
    end
    
    % KEY COMMAND: call sphereinterp_grid.m to get basis functions
    [spline_tot] = sphereinterp_grid(dlon,dlat,ax0,qparm);
    ndata = length(dlon);
    ngrid = length(spline_tot);
    
    % KEY COMMAND: call sphereinterp_est.m to perform least-squares estimation
    [dest,dest_plot,destdph_plot,destdth_plot,lam0,dlon_plot,dlat_plot,na,nb] = ...
        sphereinterp_est(spline_tot,dlon,dlat,d,dsig,ax0,rparm,pparm);
    
    if dopt==8
        %  d: real part of Z, where Z = z^2 = A + i*B
        % d2: imaginary part of Z, where Z = z^2 = A + i*B
        % (see get_1D_dataset_carl.m)
        pparm{2} = 'splitting [B]';
        [dest2,dest_plot2,destdph_plot2,destdth_plot2,lam02,dlon_plot2,dlat_plot2,na2,nb2] = ...
            sphereinterp_est(spline_tot,dlon,dlat,d2,dsig,ax0,rparm,pparm);
        
        % reconstruct the splitting vectors
        [az_deg,dt,Z] = dest2split(dest,dest2);
        figure; quiver(dlon,dlat,real(sqrt(Z)),imag(sqrt(Z)));
        
        if iwrite==1
            % estimated field at input data points
            odir = [dir_repos '/manuscripts/crichards/papers/aksplit/data/'];
            ofile = [odir stout '_est.txt'];
            fid = fopen(ofile,'w');
            for ii=1:length(dlon)
                fprintf(fid,'%12.4f%12.4f%8.1f%10.5f\n',dlon(ii),dlat(ii),az_deg(ii),dt(ii));
            end
            fclose(fid);
            
            % estimated field at regular gridpoints
            [az_deg,dt,Z] = dest2split(dest_plot,dest_plot2);
            figure; quiver(dlon_plot,dlat_plot,real(sqrt(Z)),imag(sqrt(Z)));
            
            %ftag = sprintf('%s_q%2.2i_q%2.2i_ir%2.2i_id%2.2i',slabel,qmin,qmax,ropt,dopt);
            ofile = [odir stout '_est_regular.txt'];
            fid = fopen(ofile,'w');
            for ii=1:length(dlon_plot)
                fprintf(fid,'%12.4f%12.4f%8.1f%10.5f\n',dlon_plot(ii),dlat_plot(ii),az_deg(ii),dt(ii));
            end
            fclose(fid);
        end
        
        error('stopping here for splitting');
      
    end
    
    disp('  ');
    disp(sprintf('Number of observations, ndata = %i',ndata));
    disp(sprintf('Number of basis functions, ngrid = %i',ngrid));
    disp('For testing purposes, try decreasing one of these:');
    disp(sprintf('  qmax = %i, the densest grid for basis functions',qmax));
    disp(sprintf('  nx = %i, the grid density for plotting',nx));
    disp(sprintf('  ndata = %i, the number of observations (or ax0)',ndata));
    
    % compute magnitude of surface gradient, then convert to a slope in degrees
    % note: d is in units of km, so the earth radius must also be in km
    th_plot = (90 - dlat_plot)*pi/180;
    destG_plot = sqrt( destdth_plot.^2 + (destdph_plot ./ sin(th_plot)).^2 );
    destGslope_plot = atan(destG_plot / 6371) * 180/pi;
    
    figure; scatter(dlon_plot,dlat_plot,4^2,destGslope_plot,'filled');
    axis(ax0); title('Slope of surface, degrees'); colorbar;
end

% optional: threshold the plotting field to eliminate unphysical values
% (Alaska offshore moho)
%dest_plot(dest_plot <= 11) = 11;

%----------------------------------------------------------------
% WRITE FILES

if and(iwavelet==1,iwrite==1)
    
    ftag = sprintf('%s_q%2.2i_q%2.2i_ir%2.2i_id%2.2i',slabel,qmin,qmax,ropt,dopt);
    %flab = [dir_output slabel '_' stqtag{1} '_' sprintf('ic%2.2i_im%2.2i',idata,sub_opt) ];
    flab = [dir_output ftag];
    disp('writing files with tag:'); disp(flab);
    
    nplot = length(dest_plot);
    
    % data and estimated field
    fid = fopen([flab '.dat'],'w');
    %stfmt = '%12.6f%12.6f%10.3f%10.3f%10.3f\n';
    stfmt = '%18.8e%18.8e%18.8e%18.8e%18.8e\n';
    for ii=1:ndata
        fprintf(fid,stfmt,dlon(ii),dlat(ii),d(ii),dest(ii),dsig(ii));
    end
    fclose(fid);

    % estimated field for a regular grid -- includes derivative fields, too
    fid = fopen([flab '_plot.dat'],'w');
    stfmt = '%18.8e%18.8e%18.8e%18.8e%18.8e%18.8e%18.8e\n';
    for ii=1:nplot
        fprintf(fid,stfmt,dlon_plot(ii),dlat_plot(ii),dest_plot(ii),...
            destdph_plot(ii),destdth_plot(ii),destG_plot(ii),destGslope_plot(ii));
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
    
    % save simple matlab file
    X = reshape(dlon_plot,na,nb);
    Y = reshape(dlat_plot,na,nb);
    Z = reshape(dest_plot,na,nb);
    figure; msize = 6^2; hold on;
    pcolor(X,Y,Z); shading flat; colorbar
    scatter(dlon,dlat,msize,d,'filled');
    scatter(dlon,dlat,msize,'ko');
    save(flab,'X','Y','Z');
    
    % optional: write data constraints
    error
    write_surf_data(flab,[],dlon,dlat,d,dsig,isource,ftags,szone);
    
    % testing output figure
    if 0==1
        % multiple values on the western SAF (vertical fault)
        dx = load('/home/carltape/MOHO/WAVELET/MATLAB_EST/maricopa_maricopa_basement_q05_q11_ir03_id03.dat');
        [~,ix] = sort(dx(:,3),'ascend');
        figure; scatter(dx(ix,1),dx(ix,2),4^2,dx(ix,3),'filled'); colorbar;
        [~,ix] = sort(dx(:,3),'descend');
        figure; scatter(dx(ix,1),dx(ix,2),4^2,dx(ix,3),'filled'); colorbar;
    end
end

%========================================================
