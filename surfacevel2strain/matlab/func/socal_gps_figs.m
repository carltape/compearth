%
% socal_gps_figs.m
%
% Plots test figures associated with Tape et al. (GJI 2009).
%
% called by surfacevel2strain_figs.m
%

disp('socal_gps_figs.m: additional figures for socal');
disp('  THESE ARE TEST FIGURES THAT WILL BE RE-PLOTTED GMT.');

% output directory for writing files
%odir = '/home/carltape/gmt/gps_figs/data_for_plots/';
odir = dir_output;

% file prefix
prefix = ['multiscale_' slabel '_d' sprintf('%2.2i',dopt) ];

% color scale for velocity
% first two rows are for overal and secular field
% remaining rows are for residual fields
cmat2 = [ repmat([-35 5],2,1) ; repmat([-5 5],nump-2,1)];

figure; nc=3; nr=nump;
mm = counter2(5,3); im = 1;

% velocity field: Ve and Vs as MULTI-SCALE
cfac = 0.8;
vpmask = ones(nplot,1);
for kk=1:2
    if kk==1, fvec = fe; stlab = 'veast'; end
    if kk==2, fvec = fs; stlab = 'vsouth'; end

    % plot spline gridpoints by order q
    for ip=1:nump
        inds = [iqr(ip,1) : iqr(ip,2)];

        if imask==1
            pinds = inds;
            if ip == 1, pinds = inds_qmax; end
            [m_ikeep,m_inum] = spline_thresh(G0_plot(:,pinds),m_ntrsh,m_qtrsh,[1 0]);
            vpmask = NaN * ones(nplot,1); vpmask(m_ikeep) = 1;
        end

        vplot = G0_plot(:,inds)*fvec(inds)*1e3 .* vpmask;
        vmax = cfac*max(abs(vplot));

        %subplot(nr,nc,ip+1); hold on;
        subplot(nr,nc,mm(im)); im=im+1; hold on;
        [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,vplot,npts,'cubic');
        pcolor(X,Y,Z); shading flat;
        %plot(glon(inds),glat(inds),'k.');
        axis equal, axis(ax1); caxis(cmat2(ip,1:2));
        colorbar('vert');
        if kk==1, title([stqs{ip} ', ' stis{ip}]); end

        if iwrite==1
            fid = fopen([odir prefix '_' stlab '_' stqtag{ip} '.dat'],'w');
            for iz=1:nplot
                fprintf(fid,'%18.8e%18.8e%18.8e\n',dlon_plot(iz),dlat_plot(iz),vplot(iz));
            end
            fclose(fid);
        end
    end
end

% uniform vector field for plotting
numx_plot = 16;
[vlon_plot,vlat_plot] = gridvec(lonmin,lonmax,numx_plot,latmin,latmax);
nvplot = length(vlon_plot);
vG0_plot = spline_vals_mat(spline_tot, vlon_plot, vlat_plot);

% plot uniform vector fields for the estimations
s = 1;
vpmask = ones(nvplot,1);
for ip=1:nump
    inds = [iqr(ip,1) : iqr(ip,2)];

    if imask==1
        pinds = inds;
        if ip == 1, pinds = inds_qmax; end
        [m_ikeep,m_inum] = spline_thresh(vG0_plot(:,pinds),m_ntrsh,m_qtrsh,[1 0]);
        vpmask = NaN * ones(nvplot,1); vpmask(m_ikeep) = 1;
    end
    ve_plot = vG0_plot(:,inds)*fe(inds)*1e3;
    vs_plot = vG0_plot(:,inds)*fs(inds)*1e3;

    subplot(nr,nc,mm(im)); im=im+1; hold on;
    quiver(vlon_plot,vlat_plot,ve_plot.*vpmask,-vs_plot.*vpmask,'b');
    axis equal, axis(ax1);

    if iwrite==1
        % write uniform velocity vector field to file
        if 1==1
            % list only those vectors inside the mask
            igood = find(~isnan(vpmask)==1);
            ngood = length(igood);
            disp([num2str(ngood) ' vectors for the regular field']);
            fid = fopen([odir prefix '_vec_' stqtag{ip} '.dat'],'w');
            for jj=1:length(igood)
               k = igood(jj);
               fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',vlon_plot(k),vlat_plot(k),ve_plot(k),-vs_plot(k));
            end
            fclose(fid);
        else
            % list all vectors
            fid = fopen([odir prefix '_vec_' stqtag{ip} '.dat'],'w');
            for jj=1:nvplot
                fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',vlon_plot(jj),vlat_plot(jj),ve_plot(jj),-vs_plot(jj));
            end
            fclose(fid);
        end
    end

end
fontsize(8), orient tall, wysiwyg

% additional files for GMT plotting

% q grid indexes for the plots
fid = fopen([odir prefix '_iqvec.dat'],'w');
for jj=1:nump
    fprintf(fid,'%10i%10i\n',iqvec(jj,1),iqvec(jj,2));
end
fclose(fid);

% color scale
fid = fopen([odir prefix '_vel_colors.dat'],'w');
for jj=1:nump
    fprintf(fid,'%18.8e%18.8e\n',cmat2(jj,1),cmat2(jj,2));
end
fclose(fid);

% bounds
if iwrite==1
    fid = fopen([odir prefix '_bounds.dat'],'w');
    fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',ax1(1),ax1(2),ax1(3),ax1(4)); 
    fclose(fid);
end

if imask==1
    % mask files for each level
    % GMT mask: estimated points where you WANT to plot data
    for ip=1:nump
        pinds = [iqr(ip,1) : iqr(ip,2)];
        if ip == 1, pinds = inds_qmax; end
        [m_ikeep,m_inum] = spline_thresh(G0_plot(:,pinds),m_ntrsh,m_qtrsh,[1 0]);
        vpmask = NaN * ones(nplot,1); vpmask(m_ikeep) = 1;
        if iwrite==1
            fid = fopen([odir prefix '_interior_pts_' stqtag{ip} '.dat'],'w');
            for ii=1:num
                if vpmask(ii) ~= 1 
                    fprintf(fid,'%18.8e%18.8e\n',dlon_plot(ii), dlat_plot(ii));
                end
            end
            fclose(fid);
        end
    end
end

%--------------------------------------------
% plot the vector field residuals

figure; nc=2; nr=ceil(nump/nc);
s = 1;

% overal estimated field
ve0 = ve*1e3; vs0 = vs*1e3;
stmag = [' |vmean| = ' num2str(sprintf('%.3f', mean(sqrt(ve0.^2 + vs0.^2)))) ' mm/yr'];

subplot(nr,nc,1);
quiver(dlon,dlat,ve0,-vs0,s); axis equal, axis(ax1);
xlabel(' Latitude'); ylabel(' Longitude');
title([' Data : ' stmag]);

for ip=2:nump
    inds = [iqr(ip,1) : iqr(ip,2)];
    ve_temp = G0(:,inds)*fe(inds)*1e3;
    vs_temp = G0(:,inds)*fs(inds)*1e3;
    ve0 = ve0 - ve_temp;
    vs0 = vs0 - vs_temp;
    
    stmag = [' |vmean| = ' num2str(sprintf('%.3f', mean(sqrt(ve0.^2 + vs0.^2)))) ' mm/yr'];
    
    subplot(nr,nc,ip);
    quiver(dlon,dlat,ve0,-vs0,s); axis equal, axis(ax1);
    xlabel(' Latitude'); ylabel(' Longitude');
    title({[' Data minus (q = 0-' num2str(iqvec(ip,2)) ')'], stmag});

    if iwrite==1
        fid = fopen([odir prefix '_data_remove_' stqtag{ip} '.dat'],'w');
        for k=1:ndata
           fprintf(fid,'%18.8e%18.8e%18.8e%18.8e\n',dlon(k),dlat(k),ve0(k),-vs0(k));
        end
        fclose(fid);
    end
    
end
fontsize(9), orient tall, wysiwyg

%--------------------------------------------
% scalar quantities from tensors

% datapoints for plotting
r    = earthr*ones(nplot,1);
th   = (90 - dlat_plot)/deg;
ph   = dlon_plot/deg;

figure; nc=3; nr=nump;
im = 1;

% hand-picked color plotting: cmin/cmax/cpwr, where max(z) = cmax * 10^cpwr
% shear, rotation, dilatation
if 1==1
    cmat = [
         0     3.0  -7    0 3.0    -7     -3.0 3.0    -7
         0     0.5  -7    0 0.5    -7     -0.5 0.5    -7
         0     0.8  -7    0 0.8    -7     -0.8 0.8    -7
         0     1.0  -7    0 1.0    -7     -1.0 1.0    -7
         0     3.0  -7    0 3.0    -7     -3.0 3.0    -7
         ];
end

% write color axis limits to file
fid = fopen([odir prefix '_spline_colors.dat'],'w');
for jj=1:length(cmat(:,1))
    fprintf(fid,'%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',...
        cmat(jj,1),cmat(jj,2),cmat(jj,3),cmat(jj,4),cmat(jj,5),cmat(jj,6),cmat(jj,7),cmat(jj,8),cmat(jj,9));
end
fclose(fid);
 
for ip=1:nump
%for ip=1:1
    
    % horizontal velocity field
    inds = [iqr(ip,1) : iqr(ip,2)];
    vph  = G0_plot(:,inds)*fe(inds);
    vth  = G0_plot(:,inds)*fs(inds);

    % spatial derivatives (note minus sign for vn --> vth)
    dvphdph = G0dph_plot(:,inds)*fe(inds);
    dvthdph = G0dph_plot(:,inds)*fs(inds);
    dvphdth = G0dth_plot(:,inds)*fe(inds);
    dvthdth = G0dth_plot(:,inds)*fs(inds);
    
    % ignore gradients with respect to vertical
    %dvrdth = 0*dvrdth;
    %dvrdph = 0*dvrdph;
    
    % velocity gradient tensor
    [D, W] = vel2Lmat([r th ph], [vr vth vph], dvrtpdr, ...
        [dvrdth dvthdth dvphdth], ...
        [dvrdph dvthdph dvphdph], [1 1]);
    
    % strain: six unique components
    Drr = D(:,1);
    Drth = D(:,2);
    Drph = D(:,3);
    Dthth = D(:,5);
    Dphph = D(:,9);
    Dthph = D(:,6);

    % rotation: three unique components
    Wrth = D(:,2);
    Wrph = D(:,3);
    Wthph = D(:,6);

    %----------------------------------------------
    % compute scalar quantities from tensors (see MMA notes)

    % (1) first invariant of strain: dilatation
    dilat = Drr + Dthth + Dphph;

    % (2) shear : sqrt[I2(devD)]
    shear = sqrt( Drth.^2 + Drph.^2 + Dthph.^2 + (1/3) * ...
        (Drr.^2 + Dthth.^2 + Dphph.^2 ...
        - Drr.*Dthth - Drr.*Dphph - Dthth.*Dphph ) );     % deviatoric
    %shear = sqrt( Drth.^2 + Drph.^2 + Dthph.^2 ...
    %    - Drr.*Dthth - Drr.*Dphph - Dthth.*Dphph);      % non-deviatoric

    % (3) rotat : sqrt[-I2(W)], W = devW
    % SEE MALVERN P. 131 FOR ROTATION VECTOR DETAILS
    rotat = sqrt( Wrth.^2 + Wrph.^2 + Wthph.^2 );

    vpmask = ones(nplot,1);
    if imask==1
        pinds = inds;
        if ip == 1, pinds = inds_qmax; end
        [m_ikeep,m_inum] = spline_thresh(G0_plot(:,pinds),m_ntrsh,m_qtrsh,[1 0]);
        vpmask = NaN * ones(nplot,1); vpmask(m_ikeep) = 1;
        
        shear = shear .* vpmask;
        rotat = rotat .* vpmask;
        dilat = dilat .* vpmask;
    end
    
    %===============
    
    for kk=1:3
        switch kk
            case 1, vplot = shear;
            case 2, vplot = rotat;
            case 3, vplot = dilat;
        end
        
        % color scale
        itemp = 3*kk-2 : 3*kk;
        ctemp = cmat(ip,itemp);
        clim  = ctemp(1:2);
        cnorm = 10^ctemp(3);
    
        %subplot(nr,nc,ip+1); hold on;
        subplot(nr,nc,im); im=im+1; hold on;
        [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,vplot/cnorm,npts,'cubic');
        pcolor(X,Y,Z); shading flat;
        %plot(glon(inds),glat(inds),'k.');
        axis equal, axis(ax1); caxis(clim);
        colorbar('vert');
    end
    
    if iwrite==1
        fid = fopen([odir prefix '_strain_' stqtag{ip} '.dat'],'w');
        for iz=1:nplot
            fprintf(fid,'%18.8e%18.8e%18.8e%18.8e%18.8e\n',...
                dlon_plot(iz),dlat_plot(iz),shear(iz),rotat(iz),dilat(iz));
        end
        fclose(fid);
    end
    
end
fontsize(8), orient tall, wysiwyg

break

%========================================================
