%
% surfacevel2strain_write.m
% Carl Tape, 25-July-2008
%
% This generates various velocity and strain fields for plotting in GMT.
%
% calls xxx
% called by surfacevel2strain.m
%

%========================================================
% WRITE DATA TO FILE FOR GMT PLOTTING

% KEY: file directory and label
% EXAMPLE: /home/carltape/gmt/plates/surface_velocities/misc/socal_d01_q03_q07_b1_3D_s1_u1
% socal -- region, diven by ropt
% d01 -- type of velocity field, given by dopt
% q03_q07 -- range of grid orders used
% b1 -- type of basis function, given by sopt
% 3D -- number of components
% s1 -- type of approximation used to compute the velocity gradient, given by sopt
% u1 -- whether the field had rotation removed, given by iunrotate

if istore==1, hlab = ['_' sdopt]; else hlab = ['_fix_' stref '_' smod]; end
nlab = [slabel hlab '_' stqtag{1} '_b' num2str(basistype) '_' num2str(ndim) 'D_s' num2str(sopt) '_u' num2str(iunrotate)];
%flab = [dir_plates 'surface_velocities/misc/' nlab];
flab = [dir_output nlab];

%---------------------------------

% for GMT plotting, set zero values to <0
%veps = 1e-4; izero = find( vmag_plot < veps); vmag_plot(izero) = -veps;

% mask files for each level
if imask==1
    for ip = 1:nump-1
        qmax = iqvec2(ip,2);
        disp([' qmax = ' num2str(qmax)]);
    end
    
    for ip = 1:nump-1
        qmax = iqvec2(ip,2);
        disp([' qmax = ' num2str(qmax)]);
        
        %Dscale_deg = q_scale_deg(qmax+1);
        %m_ikeep = prepare_mask_old(dlon_plot,dlat_plot,dlon,dlat,Dfrac*Dscale_deg,Nmin,Azmin);
        
        % KEY COMMAND
        inds = [iqr2(ip,1) : iqr2(ip,2)];
        %inds = [iqr(ip+1,1) : iqr(ip+1,2)];
        [m_ikeep,sigval,sigval_cum] = prepare_mask(dlon_plot,dlat_plot,...
            G_plot(1:nplot,inds),Cm_s(inds,inds),Cm_e(inds,inds),Pcum);
        
        figure; nr=3; nc=1; sedges = [-10:0.5:10];
        subplot(nr,nc,1); N = plot_histo(sigval,sedges);
        subplot(nr,nc,2); hold on; plot( sort(sigval),'.');
        ylabel('sorted sigma values');
        subplot(nr,nc,3); hold on; plot( sigval_cum,'.');
        ylabel('cumulative sigma values');
        plot([1 nplot],Pcum*[1 1],'r--');
        orient tall, wysiwyg, fontsize(10)

        figure; hold on;
        plot(dlon_plot,dlat_plot,'bx');
        plot(dlon_plot(m_ikeep),dlat_plot(m_ikeep),'r.');
        plot(dlon,dlat,'ko');
        %plot(glon(inds_qmax),glat(inds_qmax),'k.');
        legend('possible plotting gridpoint','kept plotting gridpoint','data','location','southwest');
        title(['MASK for q = ' num2str(qmax) ', keeping ' num2str(length(m_ikeep)) ' points']);
        
        pmask = NaN * ones(nplot,1);
        pmask(m_ikeep) = 1;
    
        % GMT mask: plotting points that are MASKED
        fid = fopen([flab '_masked_pts_' num2str(ip) '.dat'],'w');
        for ii=1:num
            if pmask(ii) ~= 1
                fprintf(fid,'%18.8e%18.8e\n',dlon_plot(ii), dlat_plot(ii));
            end
        end
        fclose(fid);
    end
else
    disp(' No mask for each scale will be constructed.');
end

% q grid indexes for the plots
fid = fopen([flab '_iqvec.dat'],'w');
for jj=1:nump
    fprintf(fid,'%10i%10i\n',iqvec(jj,1),iqvec(jj,2));
end
fclose(fid);
fid = fopen([flab '_iqvec2.dat'],'w');
for jj=1:nump-1
    fprintf(fid,'%10i%10i\n',iqvec2(jj,1),iqvec2(jj,2));
end
fclose(fid);
fid = fopen([flab '_qvec.dat'],'w');
for jj=1:numq
    fprintf(fid,'%10i\n',qvec(jj));
end
fclose(fid);

% color scale for multiscale VELOCITY
% first two rows are for overall and secular field
% remaining rows are for residual fields

%cmax0 = [-2 1 -15 15 -15 15]; ctick0 = [1 5 5];
cmax0 = [-5 5 -15 15 -15 15]; ctick0 = [2 5 5];
cmax1 = [5*[-1 1] 5*[-1 1] 5*[-1 1]]; ctick1 = [2 2 2];
switch floor(dopt/10)
    case 2
    case 3
    case 5
    case 6
    case 7
    case 8
        cmax0 = [-60 60 -30 30 -30 30]; ctick0 = [30 15 15]; 
        cmax1 = [5*[-1 1] 5*[-1 1] 5*[-1 1]]; ctick1 = [5 5 5];
end
vu_col = [ repmat(cmax0(1:2),2,1) ; repmat(cmax1(1:2),nump-2,1)];
vs_col = [ repmat(cmax0(3:4),2,1) ; repmat(cmax1(3:4),nump-2,1)];
ve_col = [ repmat(cmax0(5:6),2,1) ; repmat(cmax1(5:6),nump-2,1)];
cmat2 = [ vu_col vs_col ve_col];
ctick = [ repmat(ctick0,2,1) ; repmat(ctick1,nump-2,1)];

fid = fopen([flab '_vel_colors.dat'],'w');
for jj=1:nump
    fprintf(fid,'%18.8e%18.8e%18.8e%18.8e%18.8e%18.8e\n',cmat2(jj,:));
end
fclose(fid);

fid = fopen([flab '_tick_colors.dat'],'w');
for jj=1:nump
    fprintf(fid,'%18.8e%18.8e%18.8e\n',ctick(jj,:));
end
fclose(fid);

% see surfacevel2strain_figs.m
fid = fopen([flab '_colors.dat'],'w');
for jj=1:length(cmat)
    fprintf(fid,'%12.6f%12.6f%12.6f\n',cmat(jj,1),cmat(jj,2),cmat(jj,3));
end
fclose(fid);

%==================================

% DATA: velocity field (mm/yr)
fid = fopen([flab '_vec_horz_dat.dat'],'w');
for ii=1:ndata
    fprintf(fid,[repmat('%18.8e',1,7) '\n'],...
        dlon(ii),dlat(ii),1e3*ve(ii),-1e3*vs(ii),1e3*se(ii),1e3*sn(ii),0);
end

% rotational field removed (mm/yr)
if iunrotate==1
    fid = fopen([flab '_vec_horz_before_unrotate.dat'],'w');
    for ii=1:ndata
        fprintf(fid,[repmat('%18.8e',1,7) '\n'],...
            dlon(ii),dlat(ii),1e3*ve0(ii),-1e3*vs0(ii),0,0,0);
    end
    
    fid = fopen([flab '_vec_horz_rotate.dat'],'w');
    for ii=1:ndata
        fprintf(fid,[repmat('%18.8e',1,7) '\n'],...
            dlon(ii),dlat(ii),1e3*ve_ref(ii),1e3*vn_ref(ii),0,0,0);
    end
    
    % write euler vector associated with rotation
    omega_rot_rad = elatlon(3)/(1e6*deg);      % radians per year
    write_euler_vector_psxy([flab '_unrotate'],elatlon(2),elatlon(1),omega_rot_rad);
end

% velocity field COMPONENTS (ve, vn) and estimated velocity field (mm/yr)
fid = fopen([flab '_vec_horz.dat'],'w');
for ii=1:ndata
    fprintf(fid,[repmat('%18.8e',1,6) '\n'],...
        dlon(ii),dlat(ii),1e3*ve(ii),-1e3*vs(ii),1e3*ve_est(ii),-1e3*vs_est(ii));
end
fclose(fid);
if ndim==3
    fid = fopen([flab '_vec_vert.dat'],'w');
    stfmt = '%18.8e%18.8e%18.8e%18.8e\n';
    for ii=1:ndata
        fprintf(fid,stfmt,dlon(ii),dlat(ii),1e3*vu(ii),1e3*vu_est(ii));
    end
    fclose(fid);
end

% compute multiscale estimated velocity field (mm/yr)
%  (1) at data points
%  (2) at plotting points
disp(' writing the multiscale estimated velocity field...');
for kk = 1:2
    switch kk
        case 1, Gx = G; dx = dlon; dy = dlat; ptag = ''; nx = ndata;
        case 2, Gx = G_plot; dx = dlon_plot; dy = dlat_plot; ptag = '_plot'; nx = nplot;
    end

    % compute estimated velocity filed (mm/yr) for plotting, MULTISCALE UP
    if ndim==3
        Vmat = zeros(nx,nump);
        for ip = 1:nump
            inds = [iqr(ip,1) : iqr(ip,2)]';
            Vmat(:,ip) = Gx(:,inds)*fu(inds)*1e3;
        end
        fid = fopen([flab '_vfield_up' ptag '.dat'],'w');
        stfmt = [repmat('%18.8e',1,nump+2) '\n'];
        for ii=1:nx
            fprintf(fid,stfmt,dx(ii),dy(ii),Vmat(ii,:));
        end
        fclose(fid);
    end

    % compute estimated velocity filed (mm/yr) for plotting, MULTISCALE SOUTH
    Vmat = zeros(nx,nump);
    for ip = 1:nump
        inds = [iqr(ip,1) : iqr(ip,2)]';
        Vmat(:,ip) = Gx(:,inds)*fs(inds)*1e3;
    end
    fid = fopen([flab '_vfield_south' ptag '.dat'],'w');
    stfmt = [repmat('%18.8e',1,nump+2) '\n'];
    for ii=1:nx
        fprintf(fid,stfmt,dx(ii),dy(ii),Vmat(ii,:));
    end
    fclose(fid);

    % compute estimated velocity filed (mm/yr) for plotting, MULTISCALE EAST
    Vmat = zeros(nx,nump);
    for ip = 1:nump
        inds = [iqr(ip,1) : iqr(ip,2)]';
        Vmat(:,ip) = Gx(:,inds)*fe(inds)*1e3;
    end
    fid = fopen([flab '_vfield_east' ptag '.dat'],'w');
    stfmt = [repmat('%18.8e',1,nump+2) '\n'];
    for ii=1:nx
        fprintf(fid,stfmt,dx(ii),dy(ii),Vmat(ii,:));
    end
    fclose(fid);
end

% multiscale strain fields at the plotting points
disp(' writing the CUMULATIVE and INCREMENTAL multiscale estimated strain fields...');
Vmat_inc = zeros(nplot,3*(nump-1));   % by scale: q=0-5, 6-6, 7-7
Vmat_cum = zeros(nplot,3*(nump-1));   % cumulative: q=0-5, 0-6, 0-7
%disp(' writing the CUMULATIVE and INCREMENTAL multiscale estimated rotation vectors...');
%Wmat_inc = zeros(nplot,3*(nump-1));   % by scale: q=0-5, 6-6, 7-7
%Wmat_cum = zeros(nplot,3*(nump-1));   % cumulative: q=0-5, 0-6, 0-7

for ip = 1:nump-1
    inds1 = [iqr(ip+1,1) : iqr(ip+1,2)]';   % incremental
    inds2 = [iqr2(ip,1) : iqr2(ip,2)]';     % cumulative

    for kk = 1:2
        if kk==1, inds = inds1; end
        if kk==2, inds = inds2; end
    
        % vertical component of v-field, and spatial derivatives
        if ndim==3
            vr0 = G_plot(:,inds)*fu(inds);
            dvrdth0 = Gdth_plot(:,inds)*fu(inds);
            dvrdph0 = Gdph_plot(:,inds)*fu(inds);
        else
            vr0      = zeros(num,1);
            dvrdth0  = zeros(num,1);
            dvrdph0  = zeros(num,1);
        end
        vth0 = G_plot(:,inds)*fs(inds);
        vph0 = G_plot(:,inds)*fe(inds);
        dvphdph0 = Gdph_plot(:,inds)*fe(inds);
        dvthdph0 = Gdph_plot(:,inds)*fs(inds);
        dvphdth0 = Gdth_plot(:,inds)*fe(inds);
        dvthdth0 = Gdth_plot(:,inds)*fs(inds);

        % KEY: compute the velocity gradient
        [D0, W0] = vel2Lmat([r th ph], [vr0 vth0 vph0], dvrtpdr, ...
            [dvrdth0 dvthdth0 dvphdth0], ...
            [dvrdph0 dvthdph0 dvphdph0], sopt);
        Drr0 = D0(:,1); Drth0 = D0(:,2); Drph0 = D0(:,3);
        Dthth0 = D0(:,5); Dthph0 = D0(:,6); Dphph0 = D0(:,9);
        Wrth0 = W0(:,2); Wrph0 = W0(:,3); Wthph0 = W0(:,6);

        % components of strain and shear
        a1 = Drth0.^2 + Drph0.^2  + Dthph0.^2 ;
        a2 = Drr0.^2  + Dthth0.^2 + Dphph0.^2 ;
        a3 = Drr0.*Dthth0 + Drr0.*Dphph0 + Dthth0.*Dphph0 ;

        % scalar quantities
        dilat0 = Drr0 + Dthth0 + Dphph0;
        rotat0 = sqrt( Wrth0.^2 + Wrph0.^2 + Wthph0.^2 );
        %strain0 = sqrt( 2*Drth0.^2 + 2*Drph0.^2 + 2*Dthph0.^2 + Drr0.^2 + Dthth0.^2 + Dphph0.^2 );
        strain0 = sqrt( 2*a1 + a2 );                     % double-dot of strain-rate
        shear0  = sqrt( 2*a1 + (2/3)*a2 - (2/3)*a3 );    % double-dot of deviatoric strain-rate
        
        % KEY: do you want "strain" or "shear"
        ipinds = (3*ip-2):(3*ip);
        if kk==1
            Vmat_inc(:,ipinds) = [ strain0 rotat0 dilat0 ];
            %Wmat_inc(:,ipinds) = [ -Wthph0  Wrph0 -Wrth0];
        else
            Vmat_cum(:,ipinds) = [ strain0 rotat0 dilat0 ];
            %Wmat_cum(:,ipinds) = [ -Wthph0  Wrph0 -Wrth0];
        end
    end
    
end
stfmt = [repmat('%18.8e',1,3*(nump-1)+2) '\n'];
fid = fopen([flab '_multi_strain_inc.dat'],'w');
for ii=1:nplot
    fprintf(fid,stfmt,dlon_plot(ii),dlat_plot(ii),Vmat_inc(ii,:));
end
fclose(fid);
fid = fopen([flab '_multi_strain_cum.dat'],'w');
for ii=1:nplot
    fprintf(fid,stfmt,dlon_plot(ii),dlat_plot(ii),Vmat_cum(ii,:));
end
fclose(fid);

% components of strain-rate tensor
fid = fopen([flab '_Dtensor_6entries.dat'],'w');
fprintf(fid,'%18s%18s%18s%18s%18s%18s\n','Drr','Drth','Drph','Dthth','Dthph','Dphph');
for ii=1:nplot
    fprintf(fid,'%18.8e%18.8e%18.8e%18.8e%18.8e%18.8e\n',...
        Drr_all(ii),Drth_all(ii),Drph_all(ii),Dthth_all(ii),Dthph_all(ii),Dphph_all(ii));
end
fclose(fid);
% components of rotation-rate tensor
fid = fopen([flab '_Wtensor_3entries.dat'],'w');
fprintf(fid,'%18s%18s%18s\n','Wrth','Wrph','Wthph');
for ii=1:nplot
    fprintf(fid,'%18.8e%18.8e%18.8e\n',...
        Wrth_all(ii),Wrph_all(ii),Wthph_all(ii));
end
fclose(fid);

% fid = fopen([flab '_multi_wvec_inc.dat'],'w');
% for ii=1:nplot
%     fprintf(fid,stfmt,dlon_plot(ii),dlat_plot(ii),Wmat_inc(ii,:));
% end
% fclose(fid);
% fid = fopen([flab '_multi_wvec_cum.dat'],'w');
% for ii=1:nplot
%     fprintf(fid,stfmt,dlon_plot(ii),dlat_plot(ii),Wmat_cum(ii,:));
% end
% fclose(fid);

% multiscale residial field at datapoints -- lenght unit is MILLIMETERS
disp(' writing the multiscale estimated RESIDUAL velocity field...');
Vmat = zeros(ndata,3*(nump-1));
for ip = 1:nump-1
    inds = [iqr2(ip,1) : iqr2(ip,2)]';
    ipinds = (3*ip-2):(3*ip);
    if ndim==3, vu_temp = (vu - G(:,inds)*fu(inds))*1e3;
    else vu_temp = zeros(ndata,1); end
    Vmat(:,ipinds) = [ vu_temp (vs - G(:,inds)*fs(inds))*1e3 (ve - G(:,inds)*fe(inds))*1e3 ];
end
fid = fopen([flab '_vfield_residual.dat'],'w');
stfmt = [repmat('%18.8e',1,3*(nump-1)+5) '\n'];
for ii=1:ndata
    fprintf(fid,stfmt,dlon(ii),dlat(ii),su(ii)*1e3,sn(ii)*1e3,se(ii)*1e3,Vmat(ii,:));
end
fclose(fid);

% write vertical field and (overall) residual to file -- in MILLIMETERS
if ndim==3
    write_gps_psxy_vert([flab],dlon,dlat,vu*1e3,su*1e3);
    write_gps_psxy_vert([flab '_residual'],dlon,dlat,vu_res*1e3,su*1e3);
end
    
%-----------------------------------

% estimated velocity MAGNITUDE (mm/yr), dilatation, strain, rotation
% ALSO components of rotation vector
fid = fopen([flab '_strain.dat'],'w');
stfmt = [repmat('%18.8e',1,15) '\n'];
for ii = 1:num
    fprintf(fid,stfmt,...
        dlon_plot(ii), dlat_plot(ii), 1e3*vmag_plot(ii), ...
        dilat_all(ii,end), strain_all(ii,end), ...
        rotat_all(ii,end), lammax_all(ii,end), ...
        1e3*vr(ii,end), 1e3*vth(ii,end), 1e3*vph(ii,end), ...
        w_local(1,ii), w_local(2,ii), w_local(3,ii),...
        ws_all(ii), wt_all(ii));
end
fclose(fid);

% GMT mask: plotting points that are MASKED
fid = fopen([flab '_masked_pts.dat'],'w');
for ii=1:num
    if pmask(ii) ~= 1
        fprintf(fid,'%18.8e%18.8e\n',dlon_plot(ii), dlat_plot(ii));
    end
end
fclose(fid);

% write euler vectors to file
%fscale = 2;
%if and(dopt >= 20, dopt < 30), fscale = 10; end
%if and(dopt >= 50, dopt < 60), fscale = 10; end  
write_euler_vector_psxy(flab,elon,elat,omega_rad);

% write bounds to file
fid = fopen([flab '_bounds.dat'],'w');
fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',ax1(1),ax1(2),ax1(3),ax1(4));
fclose(fid);

% spline gridpoints
fid = fopen([flab '_gridpoints.dat'],'w');
for ii=1:ngrid
    fprintf(fid,'%18.8e%18.8e\n',glon(ii),glat(ii));
end
fclose(fid);

% spline gridpoints by scale
for iq=1:numq
    q = qvec(iq);
    fid = fopen([flab '_gridpoints_q' sprintf('%2.2i',q) '.dat'],'w');
    imin = ifr(iq,1);
    imax = ifr(iq,2);
    for ii=imin:imax
        fprintf(fid,'%18.8e%18.8e%6i\n',spline_tot(ii,1),spline_tot(ii,2),spline_tot(ii,3));
    end
    fclose(fid);
end

% max q for each plotting point
fid = fopen([flab '_qmax_plotting.dat'],'w');
for ii=1:nplot
    fprintf(fid,'%18.8e%18.8e%6i\n',dlon_plot(ii),dlat_plot(ii),qmax_scale(ii));
end
fclose(fid);

% for rotation fields, write out some euler poles
if any(floor(dopt/10) == [6 7])    % microplate
    pdata = [eb1m_lon eb1m_lat ; eb1am_lon eb1am_lat ; eb2m_lon eb2m_lat ; eb2am_lon eb2am_lat ];
    write_xy_points([flab '_epole_mean'],pdata(:,1),pdata(:,2));
end
if any(floor(dopt/10) == 2)         % pure rotation
    pdata = [eb1m_lon eb1m_lat ; eb1am_lon eb1am_lat ];
    write_xy_points([flab '_epole_mean'],pdata(:,1),pdata(:,2));
end

%========================================================
