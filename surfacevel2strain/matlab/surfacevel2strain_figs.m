%
% surfacevel2strain_figs.m
%
% Script for plotting various figures.
%
% calls xxx
% called by surfacevel2strain.m
%

if ifigs1==1
    disp('surfacevel2strain_figs.m: plotting with ifigs1==1');
    
    % velocity vector field: data, fit, and residuals
    figure; nr=2; nc=1;

    if ndata <= 4000
        subplot(nr,nc,1); hold on;
        quiver(dlon,dlat,ve,-vs,1,'b');
        quiver(dlon,dlat,ve_est,-vs_est,1,'r');
        if iplot_fault==1, plot(lonseg,latseg,'r','linewidth',2); end
        legend(' data field', ' estimated field');
        axis equal, axis(ax0);

        subplot(nr,nc,2); hold on;
        quiver(dlon,dlat,ve_res,-vs_res,'b');
        if iplot_fault==1, plot(lonseg,latseg,'r','linewidth',2); end
        axis equal, axis(ax0);
        title(' residual');
        
        orient tall, wysiwyg, fontsize(9)
    end

    figure; nr=3; nc=1;
    
    %[X,Y,Z] = griddataXB(dlon,dlat,vmag*1e3,npts,'nearest');
    subplot(nr,nc,1); hold on;
    %pcolor(X,Y,Z); shading flat;
    scatter(dlon,dlat,msize,vmag*1e3,'filled');
    if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
    if iplot_fault==1, plot(lonseg,latseg,'k'); end
    axis equal, axis(ax0); %caxis(clims);
    colorbar; %colorbar('horiz');
    title(' |velocity field| (data)  (mm/yr)');
    cax = caxis;

    %[X,Y,Z] = griddataXB(dlon,dlat,vmag_est*1e3,npts,'cubic');
    subplot(nr,nc,2); hold on;
    %pcolor(X,Y,Z); shading interp;
    scatter(dlon,dlat,msize,vmag_est*1e3,'filled');
    if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
    if iplot_fault==1, plot(lonseg,latseg,'k'); end
    axis equal, axis(ax0); caxis(cax);
    colorbar; %colorbar('horiz');
    title(' estimated |velocity field|  (mm/yr)');

    %[X,Y,Z] = griddataXB(dlon,dlat,vmag_res*1e3,npts,'cubic');
    subplot(nr,nc,3); hold on;
    %pcolor(X,Y,Z); shading interp;
    scatter(dlon,dlat,msize,vmag_res*1e3,'filled');
    if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
    if iplot_fault==1, plot(lonseg,latseg,'k'); end
    axis equal, axis(ax0); caxis([-1 1]*0.5*max(abs(vmag_res*1e3)));
    colorbar; %colorbar('horiz');
    title(' residual |velocity field|  (mm/yr)');
    
    orient tall, wysiwyg, fontsize(9)

    if ispheroidal==1
        disp(' plotting W, U, and V at different scales ...');
        cfac = 0.8;
        for kk = abs(ndim-4):3
            if kk==1, fvec = fU; end
            if kk==2, fvec = fV; end
            if kk==3, fvec = fW; end
            figure; nc=2; nr=ceil(nump/nc);
            for ip=1:nump
                inds = [iqr(ip,1) : iqr(ip,2)];
                cplot = G(:,inds)*fvec(inds)*1e3;
                vmax = cfac*max(abs(cplot));

                subplot(nr,nc,ip); hold on;
                %[X,Y,Z] = griddataXB(dlon,dlat,cplot,npts,'cubic');
                %pcolor(X,Y,Z); shading flat;
                scatter(dlon,dlat,msize,cplot,'filled');
                if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
                if iplot_fault==1, plot(lonseg,latseg,'k'); end
                axis equal, axis(ax0);
                %caxis([-1 1]*vmax);
                colorbar('vert');
                title([stks2{kk} ': ' stqs{ip} ', ' stis{ip}]);
            end
            orient tall, wysiwyg, fontsize(8)
        end
        
        disp(' plotting norm_vT and norm_vS (LATER: at different scales) ...');
        cfac = 0.8;
        for kk = 1:2
            figure; hold on;
            if kk==1
                [X,Y,Z] = griddataXB(dlon,dlat,norm_vS,npts,'cubic');
                pcolor(X,Y,Z); shading flat;
                %scatter(dlon, dlat, msize, norm_vS,'filled');
                if iplot_fault==1, plot(lonseg,latseg,'k'); end
                quiver( dlon, dlat, vSe, -vSs, 'k');
            else
                [X,Y,Z] = griddataXB(dlon,dlat,norm_vT,npts,'cubic');
                pcolor(X,Y,Z); shading flat;
                %scatter(dlon, dlat, msize, norm_vT,'filled');
                if iplot_fault==1, plot(lonseg,latseg,'k'); end
                quiver( dlon, dlat, vTe, -vTs, 'k');
            end

            axis equal, axis(ax0); colorbar('vert');
            title([stks3{kk} ': ' stqs{1} ', ' stis{1}]);
            orient tall, wysiwyg, fontsize(8)
        end
        
        % velocity field: ve and vs MULTI-SCALE (including the overall field)
        cfac1 = 0.8;
        figure(80); nc=2; nr=ceil(nump/nc);
        figure(81); nc=2; nr=ceil(nump/nc);
        for ip=1:nump
            inds0 = [iqr(ip,1) : iqr(ip,2)];
            inds = [inds0 inds0+ngrid];
            cplot1 = Gmode(1:ndata,inds)*fmode(inds)*1e3;
            cplot2 = Gmode(ndata+1:2*ndata,inds)*fmode(inds)*1e3;
            %vmax = cfac*max(abs(cplot));

            figure(80); subplot(nr,nc,ip); hold on;
            %[X,Y,Z] = griddataXB(dlon,dlat,cplot1,npts,'cubic');
            %pcolor(X,Y,Z); shading flat;
            scatter(dlon,dlat,msize,cplot1,'filled');
            if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
            if iplot_fault==1, plot(lonseg,latseg,'k'); end
            axis equal, axis(ax0); colorbar('vert');
            title([stks1{kk} ': ' stqs{ip} ', ' stis{ip}]);
            
            figure(81); subplot(nr,nc,ip); hold on;
            %[X,Y,Z] = griddataXB(dlon,dlat,cplot2,npts,'cubic');
            %pcolor(X,Y,Z); shading flat;
            scatter(dlon,dlat,msize,cplot2,'filled');
            if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
            if iplot_fault==1, plot(lonseg,latseg,'k'); end
            axis equal, axis(ax0); colorbar('vert');
            title([stks1{kk} ': ' stqs{ip} ', ' stis{ip}]);
        end
        figure(80), orient tall, wysiwyg, fontsize(8)
        figure(81), orient tall, wysiwyg, fontsize(8)
        
    else
        % velocity field: ve and vs MULTI-SCALE (including the overall field)
        cfac = 0.8;
        for kk = kmin:kmax
            if kk==1, fvec = fu; end
            if kk==2, fvec = fs; end
            if kk==3, fvec = fe; end
            figure; nc=2; nr=ceil(nump/nc);
            for ip=1:nump
                inds = [iqr(ip,1) : iqr(ip,2)];
                cplot = G(:,inds)*fvec(inds)*1e3;
                vmax = cfac*max(abs(cplot));

                subplot(nr,nc,ip); hold on;
                %[X,Y,Z] = griddataXB(dlon,dlat,cplot,npts,'cubic');
                %pcolor(X,Y,Z); shading flat;
                scatter(dlon,dlat,msize,cplot,'filled');
                if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
                if iplot_fault==1, plot(lonseg,latseg,'k'); end
                axis equal, axis(ax0);
                %caxis([-1 1]*vmax);
                colorbar('vert');
                title([stks1{kk} ': ' stqs{ip} ', ' stis{ip}]);
            end
            orient tall, wysiwyg, fontsize(8)
        end
    end
    
end  % ifigs1

if ispheroidal==1
    disp('surfacevel2strain_figs.m: spheroidal-toroidal decomposition');
    disp('   has not yet been extended to the regular plotting mesh.');
    disp('   EXIT HERE');
    break
end

%========================================================
% UNIFORM MESH FOR PLOTTING SCALAR FIELDS (from irregular data) --
% MASK OUT AREAS OF POOR DATA COVERAGE
% (compute new design matrices with more rows but same number of columns)

% design matrix for UNIFORM plotting grid of scalar quantities
%numx_plot = 60;
numx_plot = npts;
[dlon_plot,dlat_plot] = gridvec(lonmin,lonmax,numx_plot,latmin,latmax);
nplot = length(dlon_plot);

% design matrix for UNIFORM plotting grid of vector quantities
%numx_plot2 = 20;
%[dlon_plot2,dlat_plot2] = gridvec(lonmin,lonmax,numx_plot2,latmin,latmax);
%nplot2 = length(dlon_plot2);

%[iplot_up_i, iplot_south_i, iplot_east_i] = subindexing(nplot,ndim,opts_index);

% KEY: base design matrix for plotting
disp('  '); disp('Constructing the base design matrix for plotting...');
if basistype == 1
    [G_plot, Gdph_plot, Gdth_plot] = dogsph_vals_mat(spline_tot, dlon_plot, dlat_plot, {3});
end
if basistype == 2
    [G_plot, Gdph_plot, Gdth_plot] = spline_vals_mat(spline_tot, dlon_plot, dlat_plot);
end

if imask == 1
    % KEY COMMAND -- smaller number means larger mask
    switch floor(dopt/10)
        case 0
            % socal default
            Pcum = 0.30; if qmax == 8, Pcum = 0.45; end
            if dopt==4
                Pcum = 0.30; if qmax == 8, Pcum = 0.20; end
            end
        case 1, Pcum = 0.35;
        case 2, Pcum = 0.15;
        case 3, Pcum = 0.60;
        case 5, Pcum = 0.60;
        case 6, Pcum = 0.80;
        case 7, Pcum = 0.60;
        case 8, Pcum = 0.50;
    end
    [m_ikeep,sigval,sigval_cum] = prepare_mask(dlon_plot,dlat_plot,G_plot,Cm_s,Cm_e,Pcum);

    figure; nr=3; nc=1; sedges = [-10:0.5:10];
    subplot(nr,nc,1); N = plot_histo(sigval,sedges);
    subplot(nr,nc,2); hold on; plot( sort(sigval),'.');
    ylabel('sorted sigma values');
    subplot(nr,nc,3); hold on; plot( sigval_cum,'.');
    plot([1 nplot],Pcum*[1 1],'r--');
    ylabel('cumulative sigma values');
    orient tall, wysiwyg, fontsize(10)
    
    figure; hold on;
    plot(dlon_plot,dlat_plot,'bx');
    plot(dlon_plot(m_ikeep),dlat_plot(m_ikeep),'r+');
    plot(dlon,dlat,'ko');
    if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
    legend('possible plotting gridpoint','kept plotting gridpoint','data','location','southwest');
    
    % initialize the mask
    pmask = NaN * ones(nplot,1);
    pmask(m_ikeep) = 1;

else
    m_ikeep = [1:nplot]';
    pmask = ones(nplot,1);
end

% indices of points that are unmasked
igood = find(isnan(pmask)==0);
fprintf('%i out of %i plotting points are unmasked\n',length(igood),nplot);

% if imask == 1
%     % compute mask to white-out points that fall outside a region of
%     % (thresholded) spline coverage; this is based on using the qmax splines
% 
%     disp('  '); disp(' threshold the uniform datapoints for the mask');
% 
%     % KEY COMMANDS -- PERHAPS CHOOSE THESE BASED ON THE SYNTHETIC FIELDS,
%     % WHERE YOU CAN MASK OUT THE EDGE EFFECTS OF THE ESTIMATED FIELDS
%     Nmin = 3; Azmin = 3; Dfrac = 1.5;
%     
%     %if qmax == 8, Nmin = 3; Azmin = 3; Dfrac = 1.50; end
%     %if qmax == 7, Nmin = 3; Azmin = 3; Dfrac = 1.25; end
%     %if qmax == 6, Nmin = 3; Azmin = 3; Dfrac = 1.0; end
%     %if qmax == 5, Nmin = 3; Azmin = 3; Dfrac = 1.0; end
%     m_ikeep = prepare_mask_old(dlon_plot,dlat_plot,dlon,dlat,Dfrac*Dscale_deg,Nmin,Azmin);
%     m_imask = setdiff([1:nplot]',m_ikeep);
% 
%     %m_qtrsh = 0.5;    % evaluations must be >= qtrsh
%     %m_ntrsh = 1;      % number of evaluations that are >= qtrsh
%     %[m_ikeep,m_inum] = spline_thresh(G_plot(:,[ifr(7,1):ifr(7,2)]),m_ntrsh,m_qtrsh,[1 0]);
%     %[m_ikeep,m_inum] = spline_thresh(G_plot(:,inds_qmax),m_ntrsh,m_qtrsh,[1 0]);
%     %[m_ikeep,m_inum] = spline_thresh_2(G_plot(:,inds_qmax),m_ntrsh,m_qtrsh,glon,glat,dlon_plot,dlat_plot,[1 1 1]);
% 
%     figure; hold on;
%     plot(dlon_plot,dlat_plot,'b+');
%     plot(dlon_plot(m_ikeep),dlat_plot(m_ikeep),'r.');
%     plot(dlon,dlat,'ko');
%     %plot(glon(inds_qmax),glat(inds_qmax),'k.');
%     legend('possible plotting gridpoint','kept plotting gridpoint','data','location','southwest');
%     %legend('possible plotting gridpoint','kept plotting gridpoint',...
%     %    'data',['spline gridpoint for qmax = '
%     %    num2str(qmax)],'location','southwest');
%     
%     if 0==1
%         % loop over thresholding options to get a feel for what is best (OBSOLETE)
%         xvec = [0:0.1:0.5];
%         yvec = [1:8];
%         nr=4; nc=2;
%         kk=1; mm1 = counter(1,nr*nc,1,10);
%         for ii=1:length(xvec)
%             for jj=1:length(yvec)
%                 m_qtrsh = xvec(ii);
%                 m_ntrsh = yvec(jj);
%                 %[m_ikeep,m_inum] = spline_thresh(G_plot(:,inds_qmax),m_ntrsh,m_qtrsh,[1 0]);
% 
%                 if mm1(kk)==1, figure; end
%                 subplot(nr,nc,mm1(kk)); kk=kk+1; hold on;
%                 plot(dlon_plot,dlat_plot,'b+');
%                 plot(dlon_plot(m_ikeep),dlat_plot(m_ikeep),'r.');
%                 plot(glon(inds_qmax),glat(inds_qmax),'k.');
%                 axis equal, axis(ax0);
%                 title(['n = ' num2str(m_ntrsh) ',  q = ' num2str(m_qtrsh)]);
%                 orient tall, wysiwyg, fontsize(8)
%             end
%         end
%         break
%     end
% 
%     % initialize the mask
%     pmask = NaN * ones(nplot,1);
%     pmask(m_ikeep) = 1;
% 
% else
%     m_ikeep = [1:nplot]';
%     pmask = ones(nplot,1);
% end

% if 0==1
%     G_plot = G;
%     Gdph_plot = Gdph;
%     Gdth_plot = Gdth;
%     dlat_plot = dlat;
%     dlon_plot = dlon;
% else
%     numx_plot = 40;
%     [dlon_plot,dlat_plot] = gridvec(lonmin,lonmax,numx_plot,latmin,latmax);
%     nplot = length(dlon_plot);
%
%     [G_plot, Gdph_plot, Gdth_plot] = spline_vals_mat(spline_tot, dlon_plot, dlat_plot, ndim);
% end
% nplot = length(dlon_plot);

%--------------------------------------------------------
% MULTISCALE ASPECTS OF THE WAVELETS USED IN THE ESTIMATION PROBLEM

qmax_scale = wavelet_min_scale(spline_tot, qtrsh_cos_dist, dlon_plot, dlat_plot);

[X,Y,Z] = griddataXB(dlon_plot,dlat_plot,qmax_scale,npts,'nearest');
figure; hold on;
pcolor(X,Y,Z); shading flat;
plot(spline_tot(:,1),spline_tot(:,2),'ko');
plot(dlon,dlat,'k.');
axis equal, axis(ax0); colorbar
title('max q wavelet seen by each plotting point');

%========================================================
%========================================================
% VELOCITY SPATIAL DERIVATIVES (STRAIN RATES)
% copied from spline_wang_D_figs.m

% USER INPUT
%sopt  = input(' Type the type of strain computation (0=full; 1=elasticity; 2=viscosity): ');
sopt = 1;

num = nplot;

% estimates of the derivative fields
% (1) assume all points have equal elevation (R = radius of sphere)
% (2) assume no velocity gradients in vertical (radial) direction
%     (WHAT IS THE RATIONALE FOR THIS?)
r       = earthr*ones(num,1);
dvrtpdr = zeros(num,3);     % unobservable: gradients w.r.t. vertical

% vertical component of v-field, and spatial derivatives
if ndim==3
    vr     = G_plot*fu;
    dvrdth = Gdth_plot*fu;
    dvrdph = Gdph_plot*fu;
else
    vr      = zeros(num,1);
    dvrdth  = zeros(num,1);
    dvrdph  = zeros(num,1);
end

% horizontal velocity field
th  = (90 - dlat_plot)/deg;
ph  = dlon_plot/deg;
vth = G_plot*fs;
vph = G_plot*fe;
vmag_plot = sqrt( vth.^2 + vph.^2 );

if 0==1
    [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,vmag_plot*1e3,npts,'nearest');
    figure; hold on;
    pcolor(X,Y,Z); shading flat;
    if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
    axis equal, axis(ax0); %caxis(clims);
    colorbar('horiz'); title(' |velocity field| (data)  (mm/yr)');
    cax = caxis;
end

% spatial derivatives
dvphdph = Gdph_plot*fe;
dvthdph = Gdph_plot*fs;
dvphdth = Gdth_plot*fe;
dvthdth = Gdth_plot*fs;

% KEY: compute the velocity gradient and the curl
[D, W] = vel2Lmat([r th ph], [vr vth vph], dvrtpdr, ...
    [dvrdth dvthdth dvphdth], ...
    [dvrdph dvthdph dvphdph], sopt);

% magnitude of the curl
%curlmag = sqrt( curlv(:,1).^2 + curlv(:,2).^2 + curlv(:,3).^2 );

% strain: six unique components
Drr_all = D(:,1);
Drth_all = D(:,2);
Drph_all = D(:,3);
Dthth_all = D(:,5);
Dthph_all = D(:,6);
Dphph_all = D(:,9);

% rotation: three unique components
Wrth_all = W(:,2);
Wrph_all = W(:,3);
Wthph_all = W(:,6);

% split rotation into components with and without vr terms
[ws_r,ws_th,ws_ph,wt_r,wt_th,wt_ph] = vel2rot([r th ph], [vr vth vph], dvrtpdr, ...
    [dvrdth dvthdth dvphdth],[dvrdph dvthdph dvphdph]);
ws_all = sqrt( ws_r.^2 + ws_th.^2 + ws_ph.^2 );
wt_all = sqrt( wt_r.^2 + wt_th.^2 + wt_ph.^2 );

% % plot the curl
% figure; nr=2; nc=1;
% 
% subplot(nr,nc,1); hold on;
% plot(curlv(:,1),'g.'); plot(curlv(:,2),'b.'); plot(curlv(:,3),'r.');
% legend('curlv-r','curlv-th','curlv-ph');
% title(' curl of the velocity field');
%     
% subplot(nr,nc,2); hold on;
% [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,curlmag,npts,'cubic');
% pcolor(X,Y,Z); shading flat;
% plot(lonsaf,latsaf,'k');
% %quiver(dlon_plot,dlat_plot,curlv(:,3),-curlv(:,2),'k');
% axis equal, axis(ax0);
% colorbar('vert');
% 
% orient tall;  wysiwyg

%----------------------------------------------

% plot-specific operations for selecting subsets of euler poles
surfacevel2strain_evec;

%----------------------------------------------
% compute scalar quantities from tensors (see Mathematica notes and Latex notes)

% (1) first invariant of strain: dilatation
dilat_all = Drr_all + Dthth_all + Dphph_all;

% (2) strain
a1 = Drth_all.^2 + Drph_all.^2  + Dthph_all.^2 ;
a2 = Drr_all.^2  + Dthth_all.^2 + Dphph_all.^2 ;
a3 = Drr_all.*Dthth_all + Drr_all.*Dphph_all + Dthth_all.*Dphph_all ;
strain_all = sqrt( 2*a1 + a2 );                     % double-dot of strain-rate
shear_all  = sqrt( 2*a1 + (2/3)*a2 - (2/3)*a3 );    % double-dot of deviatoric strain-rate

% OLD
%strain_all = sqrt( abs( Drth_all.^2 + Drph_all.^2 + Dthph_all.^2 + (1/3) * ...
%    (Drr_all.^2 + Dthth_all.^2 + Dphph_all.^2 ...
%    - Drr_all.*Dthth_all - Drr_all.*Dphph_all - Dthth_all.*Dphph_all ) ));     % deviatoric
%strain_all = sqrt( abs( Drth_all.^2 + Drph_all.^2 + Dthph_all.^2 ...
%    - Drr_all.*Dthth_all - Drr_all.*Dphph_all - Dthth_all.*Dphph_all));      % non-deviatoric

% (3) rotat : sqrt[-I2(W)], W is already deviatoric
% --> SEE MALVERN P. 131 FOR ROTATION VECTOR DETAILS
rotat_all = sqrt( Wrth_all.^2 + Wrph_all.^2 + Wthph_all.^2 );
%rotat_all = Wthph_all;     % 2D version



%----------------------------------------------

% compute max abs eigenvalue of strain rate tensor
% D must have all REAL entries
%if and(~exist('gotmaxlam'), all(~isnan(D(:))) )
if all(~isnan(D(:)))
    disp(' computing the max eigenvalue at each point...');
    lammax_all  = zeros(nplot,1);
    lamdiff_all = zeros(nplot,1);
    lamsum_all  = zeros(nplot,1);

    D3 = zeros(3,3,nplot);
    D3(1,1,:) = Drr_all;
    D3(2,2,:) = Dthth_all;
    D3(3,3,:) = Dphph_all;
    D3(1,2,:) = Drth_all;
    D3(1,3,:) = Drph_all;
    D3(2,3,:) = Dthph_all;
    D3(2,1,:) = D3(1,2,:);
    D3(3,1,:) = D3(1,3,:);
    D3(3,2,:) = D3(2,3,:);

    for ii=1:nplot
        %disp([' ii = ' num2str(ii) ' / ' num2str(nplot)]);
        ltemp = eig(D3(:,:,ii));    % eigs or eig?
        lsort = sort(ltemp,'descend');
        
        lammax_all(ii)  = max(abs(lsort));
        lamdiff_all(ii) = abs( lsort(1) - lsort(3) );
        lamsum_all(ii)  = sum(lsort);
    end
    %gotlammax = 1;
end

%========================================================
% FIGURES

%if ifigs_socal==1
%    socal_gps_figs;
%end

%----------------------------------------------------------
% USER PARAMETERS FOR COLOR SCALES
% NOTE: This should be included within a single block of user parameters.
%
% hand-picked color plotting: cmin/cmax/cpwr, where max(z) = cmax * 10^cpwr
%
% THIS IS NEEDED FOR THE GMT PLOTS
% ==> dilatation, strain, rotation, lammax
%
%cmat = [-1 1 -6 ; 0 7 -7 ; -5 5 -7];  % taiwan, q=2-6
%cmat = [-3 3 -7 ; 0 3 -7 ; -3 3 -7];  % socal, q=2-6
%cmat = [-3 3 -8 ; 0 5 -8 ; -6 6 -8];  % socal, q=1-4
%cmat = [-3 3 -7 ; 0 4 -7 ; -4 4 -7];  % socal, q=4-11
%cmat = [-8 8 -8 ; 0 6 -8 ; -5 5 -8];  % tibet, q=2-5
%cmat = [-6 6 -7 ; 0 3 -7 ; 0 3 -7 ; 0 3 -7];  % socal, q=0-8
%cmat = [-1 1 -7 ; 0 2 -7 ; -2 2 -7];  % socal planar, q=0-8

% socal default values
%cmat = [-1.8 1.8 -7 ; 0 3.6 -7 ; 0 3.6 -7 ; 0 3.6 -7];
cmat = [-1.5 1.5 -7 ; 0 3 -7 ; 0 3 -7 ; 0 3 -7];  % romania colorscale

switch floor(dopt/10)
    case 0
        if qmax > 8, cmat = [-7.5 7.5 -7 ; 0 15 -7 ; 0 15 -7 ; 0 15 -7]; end
        if qmax == 8, cmat = [-1.8 1.8 -7 ; 0 3.6 -7 ; 0 3.6 -7 ; 0 3.6 -7]; end
        if dopt==4      % japan
            cmat = [-1 1 -7 ; 0 2 -7 ; 0 2 -7 ; 0 2 -7]; 
        end
    case 1
        if ropt==10, cmat = [-5 5 -9 ; 0 10 -9 ; 0 10 -9 ; 0 10 -9]; end
    case 2, cmat = [-1 1 -8 ; 0 2 -8 ; 0 2 -8 ; 0 2 -8];
    case 3, cmat = [-0.5 0.5 -5 ; 0 0.1 -4 ; 0 0.1 -4 ; 0 0.1 -4];
    case 5, cmat = [-0.9 0.9 -4 ; 0 0.9 -4 ; 0 0.9 -4 ; 0 0.9 -4];  
    case 6, cmat = [-1 1 -7 ; 0 2 -7 ; 0 0.2 -7 ; 0 1.4 -7];
        %cmat = [-0.7 0.7 -7 ; 0 1.4 -7 ; 0 0.2 -7 ; 0 1.4 -7];
    case 7, cmat = [-0.3 0.3 -7 ; 0 0.5 -7 ; 0 0.2 -7 ; 0 0.5 -7];
    case 8, cmat = [-10 10 -7 ; 0 20 -7 ; 0 20 -7 ; 0 20 -7];   
end
% check: disp([ [1:30]' floor([1:30]'/10) ])

% west US, wavelets
%if and(ropt == 1,qmax==6), cmat = [-1.4 1.4 -7 ; 0 1.4 -7 ; 0 1.4 -7 ; 0 1.4 -7]; end
%if and(ropt == 1,qmax==7), cmat = [-1.5 1.5 -7 ; 0 1.5 -7 ; 0 1.5 -7 ; 0  1.5 -7]; end

% asia
if dopt == 3, cmat = [-5 5 -8 ; 0 5 -8 ; 0 5 -8 ; 0 5 -8]; end

%----------------------------------------------------------

if ifigs2==1
    disp('surfacevel2strain_figs.m: plotting with ifigs2==1');
    if ispheroidal==1
        
        disp(' plotting velocities at different scales ...');
        cfac = 0.8;
        for kk=kmin:kmax
            if kk==1, fvec = fu; end
            if kk==2, fvec = fs; end
            if kk==3, fvec = fe; end
            figure; nc=2; nr=ceil(nump/nc);
            for ip=1:nump
                inds = [iqr(ip,1) : iqr(ip,2)];
                cplot = G_plot(:,inds)*fvec(inds)*1e3;
                vmax = cfac*max(abs(cplot));
                subplot(nr,nc,ip); hold on;
                [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,cplot,npts,'cubic');
                pcolor(X,Y,Z); shading flat;
                if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
                if iplot_fault==1, plot(lonseg,latseg,'k'); end
                axis equal, axis(ax0);
                %caxis([-1 1]*vmax);
                colorbar('vert');
                title([stks1{kk} ': ' stqs{ip} ', ' stis{ip}]);
            end
            orient tall, wysiwyg, fontsize(8)
        end
        
    else
        disp(' plotting velocities at different scales ...');
        cfac = 0.8;
        for kk=kmin:kmax
            if kk==1, fvec = fu; end
            if kk==2, fvec = fs; end
            if kk==3, fvec = fe; end
            figure; nc=2; nr=ceil(nump/nc);
            for ip=1:nump
                inds = [iqr(ip,1) : iqr(ip,2)];
                cplot = G_plot(:,inds)*fvec(inds)*1e3;
                vmax = cfac*max(abs(cplot));
                subplot(nr,nc,ip); hold on;
                [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,cplot,npts,'cubic');
                pcolor(X,Y,Z); shading flat;
                if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
                if iplot_fault==1, plot(lonseg,latseg,'k'); end
                axis equal, axis(ax0);
                %caxis([-1 1]*vmax);
                colorbar('vert');
                title([stks1{kk} ': ' stqs{ip} ', ' stis{ip}]);
            end
            orient tall, wysiwyg, fontsize(8)
        end

        figure; nc=2; nr=ceil(nump/nc);
        for ip=1:nump
            inds = [iqr(ip,1) : iqr(ip,2)];
            cplot_south = G_plot(:,inds)*fs(inds)*1e3;
            cplot_east = G_plot(:,inds)*fe(inds)*1e3;
            cplot = sqrt(cplot_south.^2 + cplot_east.^2);
            vmax = cfac*max(abs(cplot));
            subplot(nr,nc,ip); hold on;
            [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,cplot,npts,'cubic');
            pcolor(X,Y,Z); shading flat;
            if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
            if iplot_fault==1, plot(lonseg,latseg,'k'); end
            axis equal, axis(ax0);
            %caxis([-1 1]*vmax);
            colorbar('vert');
            title(['sqrt(ve^2 + vs^2): ' stqs{ip} ', ' stis{ip}]);
        end
        orient tall, wysiwyg, fontsize(8)
    end
    
    %---------------------------
    
    %icut = input(' Type 1 for making cross-section cuts of scalar maps, 0 otherwise: ');
    icut = 0;
    
    if icut==1
        
        % compute cuts
        % fstart_cuts: position along the discretized fault
        % range_cuts: length of cross-section, km
        fstart_cuts = [0.5];
        range_cuts = [200];
        if and(dopt >= 30, dopt < 40);
            fstart_cuts = [0.2 0.5 0.6];
            range_cuts = [100 200 300];
            %lonseg = lon_gc;
            %latseg = lat_gc;
        else
            fstart_cuts = [0.2 0.5 0.6];
            range_cuts = [100 200 300];
            %lonseg = lonsaf;
            %latseg = latsaf;
        end
        numcuts = length(fstart_cuts);
        half_range_cuts_deg = km2deg(range_cuts/2);
        
        snum = length(lonseg);
        figure; hold on;
        plot(lonseg,latseg,'k.');
        plot(lonseg(1),latseg(1),'rp','markersize',16);
        plot(ax0([1 2 2 1 1]),ax0([3 3 4 4 3]),'r');
        
        % normalize the discretized segment into a length = 1 curve
        [dvec_seg,azvec_seg] = lonlat2distaz(lonseg,latseg);
        dlen = sum(dvec_seg);
        dvec_seg_cumfrac = cumsum([0 ; dvec_seg ]) / dlen;
        
        nptcut_half = 250;
        nptcut = 2*nptcut_half - 1;
        cuts_lon = zeros(nptcut,numcuts);
        cuts_lat = zeros(nptcut,numcuts);
        cuts_dist = zeros(nptcut,numcuts);
        for ii = 1:numcuts
            [junk, icut] = min(abs(dvec_seg_cumfrac - fstart_cuts(ii)));
            loncen = lonseg(icut);
            latcen = latseg(icut);
            az = azvec_seg(ii);
            rng = half_range_cuts_deg(ii);
            [lat_right,lon_right] = track1(latcen,loncen,az+90,rng,[],'degrees',nptcut_half);
            [lat_left,lon_left] = track1(latcen,loncen,az-90,rng,[],'degrees',nptcut_half);
            cuts_lon(:,ii) = [flipud(lon_left(2:end)) ; lon_right ];
            cuts_lat(:,ii) = [flipud(lat_left(2:end)) ; lat_right ];
            
            % x-axis will be in km from the fault
            dvec_right = deg2km(cumsum([0 ; lonlat2distaz(lon_right,lat_right)]));
            dvec_left = deg2km(cumsum([0 ; lonlat2distaz(lon_left,lat_left)]));
            cuts_dist(:,ii) = [-flipud(dvec_left(2:end)) ; dvec_right ];
            
%             loncut = (-116 + 120)/(numcuts-1)*(i-1) - 120;
%             tmpcut = find(lonseg >= loncut);
%             loncut = lonseg(tmpcut(1));
%             latcut = latseg(tmpcut(1));
%             tgt_lon = (lonseg(tmpcut(1)+1) - loncut);
%             tgt_lat = (latseg(tmpcut(1)+1) - latcut);
%             tgt_norm = sqrt(tgt_lon^2 + tgt_lat^2);
%             tgt_lon = tgt_lon / tgt_norm;
%             tgt_lat = tgt_lat / tgt_norm;
%             rangecut = linspace(-2,2,80);
%             rangecut_lon = loncut - rangecut.*tgt_lat;
%             rangecut_lat = latcut + rangecut.*tgt_lon;
        end
        
        plot(cuts_lon,cuts_lat,'k'); axis equal
        
        %==============================
        
        % plot cuts at different scales
        figure; nc=2; nr=ceil((numcuts+1)/nc);
        cplot_south = G_plot(:,inds)*fs(inds)*1e3;
        cplot_east = G_plot(:,inds)*fe(inds)*1e3;
        cplot = sqrt(cplot_south.^2 + cplot_east.^2);
        vmax = cfac*max(abs(cplot));
        [X,Y,Z] = griddataXB(dlon_plot,dlat_plot,cplot,npts,'cubic');
        
        subplot(nr,nc,1); hold on;
        pcolor(X,Y,Z); shading flat;
        if iplot_qcen==1, plot(glon(inds_qmax),glat(inds_qmax),'k.'); end
        plot(lonseg,latseg,'k');
        plot(cuts_lon,cuts_lat,'k');
        axis equal, axis(ax0);
        %caxis([-1 1]*vmax);
        colorbar('vert');
        title(['sqrt(ve^2 + vs^2): ' stqs{1} ', ' stis{1}]);
        
        for ii = 1:numcuts
            rangecut_lon = cuts_lon(:,ii);
            rangecut_lat = cuts_lat(:,ii);
            rangecut_dist = cuts_dist(:,ii);

            if basistype == 1       % spherical wavelets
                [G_cut, Gdph_cut, Gdth_cut] = dogsph_vals_mat(spline_tot, rangecut_lon, rangecut_lat, {3});
            else                    % spherical splines
                [G_cut, Gdph_cut, Gdth_cut] = spline_vals_mat(spline_tot, rangecut_lon, rangecut_lat);
            end

            subplot(nr,nc,ii+1); hold on;
            for ip = 1:nump
                inds = [iqr(ip,1) : iqr(ip,2)];
                vcut_south = G_cut(:,inds)*fs(inds)*1e3;
                vcut_east  = G_cut(:,inds)*fe(inds)*1e3;
                vcut = sqrt(vcut_south.^2 + vcut_east.^2);
                if ip == 1
                    plot(rangecut_dist, -ones(size(rangecut_dist))*10*(ip-1),'k--');
                    plot(rangecut_dist,vcut - 10*(ip-1),'r');
                else
                    plot(rangecut_dist, -ones(size(rangecut_dist))*10*(ip-2),'k--');
                    plot(rangecut_dist,vcut - 10*(ip-2));
                end
                xlabel(' Distance from fault, km');
                axis tight;
            end
        end
        orient tall, wysiwyg, fontsize(8)
        
    end  % icut
    
    disp(' plotting strain maps...');
    % model estimates, and spherical derivatives
    figure; nr=3; nc=2;
    for ik = 1:nr*nc
        switch ik
            case 1, temp = vph(:,end); stit = 'v_\phi  (m / yr)';
            case 2, temp = vth(:,end); stit = 'v_\theta  (m / yr)';
            case 3, temp = dvphdph(:,end); stit = 'd v_\phi / d \phi  (m / yr / rad)';
            case 4, temp = dvthdph(:,end); stit = 'd v_\theta / d \phi  (m / yr / rad)';
            case 5, temp = dvphdth(:,end); stit = 'd v_\phi / d \theta  (m / yr / rad)';
            case 6, temp = dvthdth(:,end); stit = 'd v_\theta / d \theta  (m / yr / rad)';
        end
        temp = temp .* pmask;
        [X,Y,Zest] = griddataXB(dlon_plot,dlat_plot,temp,npts,stype);
        cmx = max(max(abs(Zest))); disp([ik cmx]);

        subplot(nr,nc,ik); hold on;
        pcolor(X,Y,Zest); shading flat;
        %if ik > 4, caxis(cmx*[-1 1]); else caxis(clim); end
        if iplot_fault==1, plot(lonseg,latseg,'k'); end
        colorbar('vert');
        axis equal; axis(ax0);
        title(stit);
    end
    orient tall, wysiwyg, fontsize(8)
    %colormap(seis);
    
    %------------------------------------
    
    disp(' plotting components of the scalar quantity STRAIN...');
    figure; nr=3; nc=2;
    for ik = 1:nr*nc
        switch ik
            case 1, temp = a1; stit = 'a1 = T12^2 + T13^3 + T23^2';
            case 2, temp = a2; stit = 'a2 = T11^2 + T22^3 + T33^2';
            case 3, temp = a3; stit = 'a3 = T11 T22 + T11 T33 + T22 T33 ';
            case 5, temp = strain_all(:,end); stit = ' strain = (2,1,0) . a';
            case 6, temp = shear_all(:,end); stit = '  shear = (2, 2/3, -2/3) . a';
        end
        if ik ~= 4
            temp = temp .* pmask;
            [X,Y,Zest] = griddataXB(dlon_plot,dlat_plot,temp,npts,stype);
            cmx = max(max(abs(Zest))); disp([ik cmx]);

            subplot(nr,nc,ik); hold on;
            pcolor(X,Y,Zest); shading flat;
            if iplot_fault==1, plot(lonseg,latseg,'k'); end
            colorbar('vert');
            axis equal; axis(ax0);
            title(stit);
        end
    end
    orient tall, wysiwyg, fontsize(8)

    %------------------------------------
    % overall estimates for 6 scalar quantities
    figure; nr=3; nc=2;
    dmin = min(min(dilat_all(:,end))); dmax = max(max(dilat_all(:,end)));
    for ik = 1:6
        switch ik
            case 1, temp = Dthth_all(:,end); stit = ' Dthth';
            case 2, temp = Dphph_all(:,end); stit = ' Dphph';
            case 3, temp = dilat_all(:,end); stit = ' dilatation';
            case 4, temp = Dthph_all(:,end); stit = ' Dthph';
            case 5, temp = rotat_all(:,end); stit = ' rotation';
            case 6, temp = strain_all(:,end); stit = ' strain';
        end
        temp = temp .* pmask;
        [X,Y,Zest] = griddataXB(dlon_plot,dlat_plot,temp,npts,stype);
        cmx = max(max(abs(Zest)));
        stcmx = num2str(sprintf('%.4e',cmx));
        disp([num2str(ik) '  ' stcmx]);

        subplot(nr,nc,ik); hold on;
        pcolor(X,Y,Zest); shading interp;
        if iplot_fault==1, plot(lonseg,latseg,'k'); end
        if ik <= 4, caxis([dmin dmax]); end
        axis equal; axis(ax0); colorbar('vert');
        title(stit);
        pause(0.5)
    end
    %fontsize(8)
    orient tall; wysiwyg % colormap(seis);

    %------------------------------------
    % overall estimates for 3 scalar quantities (see cmat above)
    figure; nr=2; nc=2;
    for ik = 1:4
        %for ik = 4:4

        % color scale
        clim = cmat(ik,1:2);
        cnorm = 10^cmat(ik,3);
        switch ik
            case 1, temp = dilat_all(:,end); stit = ' dilatation';
            case 2, temp = strain_all(:,end); stit = ' strain';
            case 3, temp = rotat_all(:,end); stit = ' rotation';
            case 4, temp = lammax_all(:,end); stit = ' lammax';
        end
        temp = temp .* pmask;

        [X,Y,Zest] = griddataXB(dlon_plot,dlat_plot,temp/cnorm,npts,stype);
        [cmx, imax] = max(abs(Zest(:)));
        zmax = Zest(imax) * cnorm;
        stcmx = num2str(sprintf('%.4e',zmax));
        disp([num2str(ik) '  ' stcmx]);

        subplot(nr,nc,ik); hold on;
        %figure; hold on;
        pcolor(X,Y,Zest); shading interp;
        if iplot_fault==1, plot(lonseg,latseg,'k'); end
        axis equal; axis(ax0); caxis(clim);
        colorbar('vert');
        title([stit '  (max = ' stcmx ')']);
    end
    orient tall, wysiwyg, colormap(seis), fontsize(8)

    % overall estimates for 3 scalar quantities -- DEFAULT color scaling
    figure; nr=3; nc=2;
    for ik = 1:nr*nc
        switch ik
            case 1, temp = dilat_all(:,end); stit = ' dilatation';
            case 2, temp = strain_all(:,end); stit = ' strain';
            case 3, temp = rotat_all(:,end); stit = ' rotation';
            case 4, temp = lammax_all(:,end); stit = ' lammax';
            case 5, temp = lamdiff_all(:,end); stit = ' lam-diff';
            case 6, temp = lamsum_all(:,end); stit = ' lam-sum';
                
        end
        temp = temp .* pmask;
        [X,Y,Zest] = griddataXB(dlon_plot,dlat_plot,temp,npts,stype);
        subplot(nr,nc,ik); hold on;
        pcolor(X,Y,Zest); shading interp;
        if iplot_fault==1, plot(lonseg,latseg,'k'); end
        axis equal; axis(ax0);
        colorbar('vert'); title(stit);
    end
    orient tall, wysiwyg, colormap(seis), fontsize(8)

    % overall estimates for 3 scalar quantities -- DEFAULT axis scale
    % NOTE: plot only the points that are UNmasked
    figure; nr=2; nc=2;
    for ik = 1:nr*nc
        switch ik
            case 1, temp = dilat_all(igood,end); stit = ' dilatation';
            case 2, temp = strain_all(igood,end); stit = ' strain';
            case 3, temp = rotat_all(igood,end); stit = ' rotation';
            case 4, temp = lammax_all(igood,end); stit = ' lammax';
        end
        subplot(nr,nc,ik); hold on;
        plot(temp,'.');
        title({[stit ' : median, mean, std'],...
            sprintf('%.3e, %.3e, %.3e',median(temp),mean(temp),std(temp))});
    end
    orient tall, wysiwyg, fontsize(10)

    % plot component entries
    figure;
    subplot(2,1,2); hold on;
    plot(-Wthph_all(:,end),'g.');
    plot(Wrph_all(:,end),'r.');
    plot(-Wrth_all(:,end),'b.');
    plot(rotat_all(:,end),'k.');
    ax = axis;
    legend('wr = -Wthph','wth = -Wrphi','wphi = Wrth','|w|');
    subplot(4,2,1); plot(-Wthph_all(:,end),'g.'); axis(ax); title('wr = -Wthph');
    subplot(4,2,2); plot(Wrph_all(:,end),'r.'); axis(ax);  title('wth = Wrph');
    subplot(4,2,3); plot(-Wrth_all(:,end),'b.'); axis(ax);  title('wph = -Wrth');
    subplot(4,2,4); plot(rotat_all(:,end),'k.'); axis(ax);  title('|w|');
    orient tall;  wysiwyg

    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on;
    plot(Drr_all(:,end),'g.'); plot(Dthth_all(:,end),'b.'); plot(Dphph_all(:,end),'r.');
    legend('Drr','Dthth','Dphph');
    title(' diagonal elements of strain tensor');
    subplot(nr,nc,2); hold on;
    plot(Dthph_all(:,end),'g.'); plot(Drth_all(:,end),'b.'); plot(Drph_all(:,end),'r.');
    legend('Dthph','Drth','Drph');
    title(' off-diagonal elements of strain tensor');
    orient tall;  wysiwyg

end  % ifigs2

disp('exiting surfacevel2strain_figs.m');

%========================================================
