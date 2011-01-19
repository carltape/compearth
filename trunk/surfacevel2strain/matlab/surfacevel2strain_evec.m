%
% surfacevel2strain_evec.m
% Carl Tape, 12-June-2009
%
% 
%
% called by surfacevel2strain_figs.m
% calls xxx
%

% invariant quantity
rotat_radyr = sqrt( Wrth_all.^2 + Wrph_all.^2 + Wthph_all.^2 );
rotat_degMyr = rotat_radyr*1e6*deg;

% convert to euler vectors: lat, lon, omega
% KEY NOTE: m_ikeep for only the points not under the mask
Pxyz = tp2xyz(th, ph, r);
w_local = [-Wthph_all Wrph_all -Wrth_all]';  % NOTE SIGN
w_global = global2local(w_local, Pxyz, 0);
e_latlon = euler_convert(w_global*1e6*deg,1)';
elat  = e_latlon(m_ikeep,1);
elon  = e_latlon(m_ikeep,2);
omega = e_latlon(m_ikeep,3);   % deg/Myr
omega_rad = omega/(1e6*deg);   % rad/yr
[elat_anti,elon_anti] = antipode(elat,elon);

% plot euler poles and antipoles
figure; nr=2; nc=1;
%ibig = find(omega_rad > 0.2e-7);
subplot(nr,nc,1); hold on; plot(elon,elat,'b.'); title('euler poles');
%plot(elon(ibig),elat(ibig),'k.'); 
subplot(nr,nc,2); hold on; plot(elon_anti,elat_anti,'r.'); title('euler anti-poles');
%plot(elon_anti(ibig),elat_anti(ibig),'k.'); 

%=====================================================================

if and(ropt==10, floor(dopt/10) == 1)
    % select two sets of euler poles
    frac = 0.1;
    omega1 = 2.75e-9;   % qet poles and anti-poles within frac of this
    omega2 = 10e-9;     % get poles and anti-poles GREATER than this
    
    iblock1 = find( and( abs(omega_rad) >= (1-frac)*omega1, abs(omega_rad) <= (1+frac)*omega1 ) );
    iblock2 = find( abs(omega_rad) > omega2 );
    ibig = union(iblock1,iblock2);
    
    % REPLACE the full set of euler vectors with the reduced set
    elat  = elat(ibig);
    elon  = elon(ibig);
    omega = omega(ibig);   % deg/Myr
    omega_rad = omega/(1e6*deg);   % rad/yr
    [elat_anti,elon_anti] = antipode(elat,elon);
    
    figure; nr=3; nc=1;
    subplot(nr,nc,1); hold on;
    plot( dlon_plot(m_ikeep(ibig)), dlat_plot(m_ikeep(ibig)), '.');
    for ii=1:length(ibig)
        jj = m_ikeep(ibig(ii));
        text( dlon_plot(jj), dlat_plot(jj), num2str(ii));
    end
    subplot(nr,nc,2); hold on; plot( elon, elat, 'b.');
    for ii=1:length(ibig), text( elon(ii), elat(ii), num2str(ii)); end
    subplot(nr,nc,3); plot( elon_anti, elat_anti, 'r.');
    for ii=1:length(ibig), text( elon_anti(ii), elat_anti(ii), num2str(ii)); end
    orient tall, wysiwyg, fontsize(10)
    
elseif any(floor(dopt/10) == [6 7])
    % compute the best-fitting euler poles for the two-microplate field
    
    % pre-picked cluster centers with spreads
    if floor(dopt/10) == 6
        elon1_tar = -119; elat1_tar = 33.5; scirc1 = 0.5; omega1 = 1.7453e-08; frac = 0.1;
        elon2_tar = 62.5; elat2_tar = -36; scirc2 = 0.5; omega2 = 8.7266e-09; frac = 0.1;
        elon3_tar = 66; elat3_tar = -32; scirc3 = 2; omega3 = 0.2e-7;
    else
        elon1_tar = -119; elat1_tar = 33.5; scirc1 = 0.5; omega1 = 1.7453e-8; frac = 0.1;
        elon2_tar = -117.5; elat2_tar = 36; scirc2 = 0.5; omega2 = 0.87266e-8; frac = 0.1;
        elon3_tar = 65; elat3_tar = -35; scirc3 = 10; omega3 = 1e-8;
    end
    
    % find all euler vectors that are within X of these rotation values
    iblock1 = find( and( omega_rad >= (1-frac)*omega1, omega_rad <= (1+frac)*omega1 ) );
    iblock2 = find( and( omega_rad >= (1-frac)*omega2, omega_rad <= (1+frac)*omega2 ) );
    iblock3 = find( omega_rad >= omega3 );
    
    % find the points in the cluster that are within a certain radius of each point
    dcirc1 = arcdist(elat1_tar, elon1_tar, elat(iblock1), elon(iblock1));
    dcirc2 = arcdist(elat2_tar, elon2_tar, elat(iblock2), elon(iblock2));
    dcirc3 = arcdist(elat3_tar, elon3_tar, elat(iblock3), elon(iblock3));
    idcirc1 = find( dcirc1 <= scirc1 );
    idcirc2 = find( dcirc2 <= scirc2 );
    idcirc3 = find( dcirc3 <= scirc3 );
    iblock1 = iblock1(idcirc1);
    iblock2 = iblock2(idcirc2);
    iblock3 = iblock3(idcirc3);
    
    figure; subplot(2,1,1); hold on;
    plot(dlon_plot(m_ikeep(iblock1)),dlat_plot(m_ikeep(iblock1)),'ro');
    plot(dlon_plot(m_ikeep(iblock2)),dlat_plot(m_ikeep(iblock2)),'bo');
    plot(dlon_plot(m_ikeep(iblock3)),dlat_plot(m_ikeep(iblock3)),'ko');
    
    subplot(2,3,4); plot(elon(iblock1),elat(iblock1),'r.');
    subplot(2,3,5); plot(elon(iblock2),elat(iblock2),'b.');
    subplot(2,3,6); plot(elon(iblock3),elat(iblock3),'k.');
    
    % compute the mean pole for each of the two blocks
    eblock1      = latlon2xyz(elat(iblock1), elon(iblock1), earthr);
    eblock1_anti = latlon2xyz(elat_anti(iblock1), elon_anti(iblock1), earthr);
    eblock2      = latlon2xyz(elat(iblock2), elon(iblock2), earthr);
    eblock2_anti = latlon2xyz(elat_anti(iblock2), elon_anti(iblock2), earthr);
    [eb1m_lat, eb1m_lon]   = xyz2latlon(sum(eblock1,2));
    [eb1am_lat, eb1am_lon] = xyz2latlon(sum(eblock1_anti,2));
    [eb2m_lat, eb2m_lon]   = xyz2latlon(sum(eblock2,2));
    [eb2am_lat, eb2am_lon] = xyz2latlon(sum(eblock2_anti,2));
    
    % REPLACE the full set of euler vectors with the reduced set
    ibig = unique([ iblock1 ; iblock2 ; iblock3]);
    elat  = elat(ibig);
    elon  = elon(ibig);
    omega = omega(ibig);   % deg/Myr
    omega_rad = omega/(1e6*deg);   % rad/yr
    [elat_anti,elon_anti] = antipode(elat,elon);
    
    fprintf('%i euler poles kept out of %i inside the mask -- %i plotting points total\n',...
        length(ibig),length(m_ikeep),nplot);
    
%     % pick two small circle regions with near-uniform rotation
%     elon_cen1 = -120; elat_cen1 = 31; scirc1 = 0.5;
%     elon_cen2 = -115; elat_cen2 = 37; scirc2 = 0.5;
%     
%     % find the plotting points that are within a certain radius of each point
%     dcirc1 = arcdist(elat_cen1, elon_cen1, dlat_plot, dlon_plot);
%     dcirc2 = arcdist(elat_cen2, elon_cen2, dlat_plot, dlon_plot);
%     dcirc3 = arcdist(elat_cen3, elon_cen3, dlat_plot, dlon_plot);
%     idcirc1 = find( dcirc1 <= scirc1 );
%     idcirc2 = find( dcirc2 <= scirc2 );
% 
%     figure; hold on;
%     plot(dlon_plot(idcirc1),dlat_plot(idcirc1),'r.');
%     plot(dlon_plot(idcirc2),dlat_plot(idcirc2),'b.');
%     plot(elon_cen1,elat_cen1,'kp','markersize',18);
%     plot(elon_cen2,elat_cen2,'kp','markersize',18);
%     
%     % find all euler vectors that are within X of these rotation values
%     frac = 0.1;
%     omega1 = mean(rotat_degMyr(idcirc1));
%     omega2 = mean(rotat_degMyr(idcirc2));
%     iblock1 = find( and( rotat_degMyr >= (1-frac)*omega1, rotat_degMyr <= (1+frac)*omega1 ) );
%     iblock2 = find( and( rotat_degMyr >= (1-frac)*omega2, rotat_degMyr <= (1+frac)*omega2 ) );
%     plot(dlon_plot(iblock1),dlat_plot(iblock1),'ro');
%     plot(dlon_plot(iblock2),dlat_plot(iblock2),'bo');
%     
%     % plot the selected euler poles
%     %elat_b1 = elat(iblock1); elon_b1 = elon(iblock1);
%     %elat_b2 = elat(iblock2); elon_b2 = elon(iblock2);
%     %plot(elon(iblock1), elat(iblock1), 'r+', elon_anti(iblock1), elat_anti(iblock1), 'r+');
%     %plot(elon(iblock2), elat(iblock2), 'b+', elon_anti(iblock2), elat_anti(iblock2), 'b+');
%     
%     % plot the selected euler poles and pick the two clusters
%     figure; plot(elon(iblock1), elat(iblock1), 'r+');
%     elon1_tar = input(' Enter approximate longitude of desired cluster: ');
%     elat1_tar = input(' Enter approximate latitude of desired cluster: ');
%     figure; plot(elon(iblock2), elat(iblock2), 'r+');
%     elon2_tar = input(' Enter approximate longitude of desired cluster: ');
%     elat2_tar = input(' Enter approximate latitude of desired cluster: ');
%     scirc0 = input(' Enter angular radius of small circle containing the cluster of poles: ');
%     
%     % find the points in the cluster that are within a certain radius of each point
%     dcirc1 = arcdist(elat1_tar, elon1_tar, elat(iblock1), elon(iblock1));
%     dcirc2 = arcdist(elat2_tar, elon2_tar, elat(iblock2), elon(iblock2));
%     idcirc1 = find( dcirc1 <= scirc0 );
%     idcirc2 = find( dcirc2 <= scirc0 );
%     disp(sprintf('%i out of %i that are less than %.1f deg from (%.2f, %.2f)',...
%         length(idcirc1),length(dcirc1),scirc0,elon1_tar,elat1_tar));
%     disp(sprintf('%i out of %i that are less than %.1f deg from (%.2f, %.2f)',...
%         length(idcirc2),length(dcirc2),scirc0,elon2_tar,elat2_tar));
%     iblock1 = iblock1(idcirc1);
%     iblock2 = iblock2(idcirc2);
%     
%     % compute the mean pole for each block
%     eblock1      = latlon2xyz(elat(iblock1), elon(iblock1), earthr);
%     eblock1_anti = latlon2xyz(elat_anti(iblock1), elon_anti(iblock1), earthr);
%     eblock2      = latlon2xyz(elat(iblock2), elon(iblock2), earthr);
%     eblock2_anti = latlon2xyz(elat_anti(iblock2), elon_anti(iblock2), earthr);
%     [eb1m_lat, eb1m_lon]   = xyz2latlon(sum(eblock1,2));
%     [eb1am_lat, eb1am_lon] = xyz2latlon(sum(eblock1_anti,2));
%     [eb2m_lat, eb2m_lon]   = xyz2latlon(sum(eblock2,2));
%     [eb2am_lat, eb2am_lon] = xyz2latlon(sum(eblock2_anti,2));
%     axis equal; axis(ax1);
%     
%     % REPLACE the full set of euler vectors with the reduced set
%     ibig = unique([ iblock1 iblock2 ]);
%     elat  = elat(ibig);
%     elon  = elon(ibig);
%     omega = omega(ibig);   % deg/Myr
%     omega_rad = omega/(1e6*deg);   % rad/yr
%     [elat_anti,elon_anti] = antipode(elat,elon);
    
else
    
    omega_min = input(' Enter minimum value of omega for euler vectors (default = 0 rad/yr): ');
    ibig = find( abs(omega_rad) >  omega_min );
    
    % REPLACE the full set of euler vectors with the reduced set
    elat  = elat(ibig);
    elon  = elon(ibig);
    omega = omega(ibig);   % deg/Myr
    omega_rad = omega/(1e6*deg);   % rad/yr
    [elat_anti,elon_anti] = antipode(elat,elon);
    
    % compute the mean pole
    eblock1      = latlon2xyz(elat, elon, earthr);
    eblock1_anti = latlon2xyz(elat_anti, elon_anti, earthr);
    [eb1m_lat, eb1m_lon]   = xyz2latlon(sum(eblock1,2));
    [eb1am_lat, eb1am_lon] = xyz2latlon(sum(eblock1_anti,2));
    
    % plot the local poles
    if 1==1
        figure; nr=3; nc=1;
        subplot(nr,nc,1); hold on;
        plot( dlon_plot(m_ikeep(ibig)), dlat_plot(m_ikeep(ibig)), '.');
        for ii=1:length(ibig)
            jj = m_ikeep(ibig(ii));
            text( dlon_plot(jj), dlat_plot(jj), num2str(ii));
        end
        subplot(nr,nc,2); hold on; plot( elon, elat, 'b.');
        for ii=1:length(ibig), text( elon(ii), elat(ii), num2str(ii)); end
        subplot(nr,nc,3); plot( elon_anti, elat_anti, 'r.');
        for ii=1:length(ibig), text( elon_anti(ii), elat_anti(ii), num2str(ii)); end
        orient tall, wysiwyg, fontsize(10)
    end
end

%break

if 0==1
    [elat_anti,elon_anti] = antipode(elat,elon);
    for jj = 1:5
        ii = randomint(1,length(m_ikeep),1);
        kk = m_ikeep(ii);
        %ii = getsubset(dlon_plot,dlat_plot,[-116.0128 -115.9878   33.2345   33.3002]);

        % plotting the antipodes
        figure; hold on;
        plot(dlon_plot.*pmask,dlat_plot.*pmask,'k+');
        plot(elon_anti,elat_anti,'.');
        plot(-116,35,'ko','markersize',15,'markerfacecolor','r');
        plot(elon_anti(ii),elat_anti(ii),'ko','markersize',10,'markerfacecolor','w');
        plot(dlon_plot(kk),dlat_plot(kk),'rp','markersize',20,'markerfacecolor','k');
    end
end

%========================================================
