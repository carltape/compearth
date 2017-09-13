%
% run_Vhat_lib.m
%
% Compute a library of Vhat'(omega) curves that show the uniform
% distribution of moment tensors for a fixed source type (i.e., fixed lune
% point or fixed eigenvalues).
%
% W. Tape and C. Tape, 2017, GJI
% Volume in moment tensor space in terms of distance
%
% calls Vgammaomega.m
%

close all
clc
clear

deg = 180/pi;
path_Vhat;
bdir = Vhatdir_write;
ddir = strcat(bdir,'data/temp/');
pdir = strcat(bdir,'figs/temp/');
xticks = [0:30:180];
ax0 = [-1/deg 181/deg -0.02 1.4];

%--------------------------------------------------------------------------
% KEY COMMANDS
bcompute_lib = false;
bwrite_Vhat = false;
bwrite_omegac = false;
bplot_all = false;
bprint = false;

atol_omegac = 1e-12;

domega = 0.1;                   % DEGREES
omega = [0:domega:180]/deg;     % RADIANS

dgamma = 0.5;                   % DEGREES
gammavec = [0:dgamma:30]/deg;   % RADIANS
%gammavec = -[0:dgamma:30]/deg; gammavec([1 end]) = [];
%--------------------------------------------------------------------------

nomega = length(omega);
ngamma = length(gammavec);
if mean(gammavec) > 0
    gtag = 'gammap'; gtag2 = '\gamma >= 0'; gtag3 = '';
else
    gtag = 'gamman'; gtag2 = '\gamma < 0'; gtag3 = '-';
end

% integration limits
mumin = 0; mumax = pi/2;
numin = 0; numax = pi;

kmin = 1; kmax = ngamma;
%kmin = ngamma; kmax = kmin;
%kmin = 15; kmax = kmin;

% delta_deg = 10;true
% for kk=kmin:kmax
%     gamma_deg = gammavec(kk)*deg;
%     lam = lune2lam(gamma_deg,delta_deg);
%     [omegacrit,isort] = lam2omegacrit(lam);
%     omegacrit(:,4) = [];    % remove redundant critical angle
%     disp(omegacrit)
% end

%break

if bcompute_lib
    
    if bwrite_Vhat
        % list of gamma values
        ofile = sprintf('%sgammavec.dat',ddir);
        fid = fopen(ofile,'w');
        for kk=kmin:kmax
            fprintf(fid,'%8.2f  %4.4i %3i\n',gammavec(kk)*deg,round(100*gammavec(kk)*deg),kk);   
        end
        fclose(fid);

        % list of omega values
        ofile = sprintf('%somegavec.dat',ddir);
        fid = fopen(ofile,'w');
        for ii=1:length(omega), fprintf(fid,'%16.8f\n',omega(ii)*deg); end
        fclose(fid);
    end
    
    % loop over gamma values
    for kk=kmin:kmax    
        disp(sprintf('k = %i [%i to %i]',kk,kmin,kmax));
        gamma = gammavec(kk);       % RADIANS
        
         % calculate and save Vhat curves as matlab files
        if bwrite_Vhat
            fname = sprintf('%sVhat_100gamma_%4.4i.mat',ddir,round(100*gamma*deg))
            if exist(fname,'file')
                disp('file already exists');
            else
                [V,tocv,Vp,Vn,tocvp,tocvn] = Vgammaomega(gamma,omega);
                if bwrite_Vhat
                    save(fname,'omega','V','tocv','Vp','Vn','tocvp','tocvn');
                end
            end
        end
        
        % critical angles for delta = 0 Vhat curves
        if bwrite_omegac
            lam = lune2lam(gamma*deg,0);
            omegacrit_deg = lam2omegacrit(lam);
            omegacritmod_deg = omegacritmod(omegacrit_deg);
            % do not calculate Vhat for omega0 or omega4
            if length(omegacritmod_deg) > 2
                [Vpcrit,Vcrit] = Vomega_point(gamma*deg,0,omegacritmod_deg(2:end-1),atol_omegac);
                Vpcrit  = [0 ; Vpcrit(:) ; 0];
                Vcrit   = [0 ; Vcrit(:) ; 1];
            else
                Vpcrit  = [0 ; 0];
                Vcrit   = [0 ; 1];
            end

            % write to ASCII file
            fname = sprintf('%somegacrit_100gamma_%4.4i.dat',ddir,round(100*gamma*deg))
            fid = fopen(fname,'w');
            for xx=1:length(omegacritmod_deg)
                fprintf(fid,'%20.12f %20.12f %20.12f\n',omegacritmod_deg(xx),Vpcrit(xx),Vcrit(xx));   
            end
            fclose(fid);
        end
    end
    
else

    % for plotting
    warning('plotting from a different directory from where the data files are being writting');
    warning('(you might want to load the pre-saved gammavec.dat)');
    
    % tmax: max time to evaluate Vhat, seconds
    %ftag = 'temp'; tmax = 180;
    %ftag0 = 'iterated_atol08_rtol14'; tmax = 10;
    %ftag0 = 'iterated_atol10_rtol18'; tmax = 60;
    ftag0 = 'iterated_atol10_rtol18';
    if min(gammavec) < 0, tmax = 180; else tmax = 60; end
    tmax = 60;
    ftag = strcat(gtag,'_',ftag0);
    
    ddir = sprintf('%s/data/%s/',bdir,ftag);
    pdir = sprintf('%s/figs/%s/',bdir,ftag);
    %wdir = sprintf('%s/data/%s_omegacrit_atol12_rtol18_domega06/',bdir,gtag);
    disp(sprintf('plots using the files from:'));
    %disp(sprintf('INPUT omegac: %s',wdir));
    disp(sprintf('INPUT VHAT: %s',ddir));
    disp(sprintf('OUTPUT PLOTS: %s',pdir));
    
    Vtime = NaN(ngamma,nomega);
    
    for kk=kmin:kmax
        % gamma
        gamma_deg = gammavec(kk)*deg;

        % check if pre-computed Matlab file exists
        ftag = sprintf('Vhat_100gamma_%4.4i',round(gamma_deg*100));
        dfile = sprintf('%s/%s.mat',ddir,ftag);
        if ~exist(dfile,'file');
            disp('file does not exist:'); dfile
            continue
        end
        
        % vertical lines at critical angles
        delta_deg = 0;
        lam = lune2lam(gamma_deg,delta_deg);
        [omegacrit,isort] = lam2omegacrit(lam);
        omegacrit(:,4) = [];    % remove redundant critical angle
        
        % max omega value
        omegamax = max(omegacrit)/deg;
        
        % load pre-computed Vhat curve
        load(dfile);
        % NaN
        V(omega > omegamax) = NaN;
        tocv(omega > omegamax) = NaN;
        tocvn(omega > omegamax) = NaN;
        tocvp(omega > omegamax) = NaN;
        % compute Vhat prime
        domega = omega(2)-omega(1);
        nimag = norm(imag(V))/norm(real(V));
        V = real(V);
        Vprime = gradient(V,omega);
        
        Vtime(kk,:) = tocv(:)';
        
        if bplot_all
            figure; nr=2; nc=1;
            for xx=1:2
                subplot(nr,nc,xx); hold on;
                if xx==1
                    [xmat,ymat] = vertlines(omegacrit/deg,ax0(3),ax0(4));
                    plot(xmat,ymat,'k--');
                    dxx = -6;
                    ytxt = ax0(4) - 0.05*(ax0(4)-ax0(3));
                    for jj=1:4, text((omegacrit(jj)+dxx)/deg,ytxt,sprintf('\\omega_%i',jj)); end

                    h1 = plot(omega,Vprime,'b-','linewidth',1);
                    h0 = plot(omega,V,'k-','linewidth',1);
                    legend([h1 h0],'Vhat''(\omega)','Vhat(\omega)','location','northwest');
                    title(sprintf('Vhat''(\\omega) for lune point (\\gamma = %.1f^\\circ, \\delta = %.1f^\\circ)',gamma_deg,delta_deg));
                    axis equal, axis(ax0);
                    %ylabel('Vhat(\omega) and Vhat''(\omega)');
                else
                    [xmat, ymat] = vertlines(omegacrit/deg,0,tmax);
                    plot(xmat,ymat,'k--');

                    h1 = plot(omega,tocv,'k.');
                    h2 = plot(omega,tocvn,'r.');
                    h3 = plot(omega,tocvp,'c.');
                    legend([h1 h2 h3],'t(Vn)+t(Vp), sec','t(Vn), sec','t(Vp), sec','location','northwest');
                    ylabel('time to evaluate Vhat(\omega), seconds');
                    title(sprintf('total time for %i \\omega values = %.1f min (n%.1f + p%.1f)',length(omega), ...
                        sum(tocv(~isnan(tocv)))/60, ...
                        sum(tocvn(~isnan(tocvn)))/60, ...
                        sum(tocvp(~isnan(tocvp)))/60 ));
                    axis([ax0(1:2) -0.1 tmax]);
                end
                set(gca,'xtick',xticks/deg,'xticklabel',numvec2cell(xticks));
                xlabel('\omega, deg');
                %title(sprintf('Vhat''(\\omega) for lune point (\\gamma = %.1f, \\delta = %.1f) [Im/Re = %.1e]',...
                %    gamma_deg,delta_deg,nimag));
            end

            if bprint
                orient tall; print(gcf,'-dpdf',sprintf('%s/%s',pdir,ftag));
                pause(1); close all;
            end
        end  % bplot
    end
    
    % exit if you are only plotting one gamma
    if kmax-kmin==0, break; end
    
    % 2D plot of gamma and omega
    Vplot = Vtime;
    [W,G] = meshgrid(omega*deg,abs(gammavec)*deg);
    figure; hold on;
    pcolor(W,G,Vplot); shading flat; axis([0 180 0 30]);
    xlabel('\omega, deg'); ylabel('\gamma, degrees');
    caxis([0 tmax]); colorbar
    title(sprintf('time to calculate Vhat_{%s\\gamma}(\\omega), seconds',gtag3)); 
    for ii=1:length(gammavec)
        gamma0 = abs(gammavec(ii))*deg;
        delta0 = 0;
        lam = lune2lam(gamma0,delta0);
        omegac = lam2omegacrit(lam);
        omegac = omegacritmod(omegac);
        omegac(1) = [];
        plot(omegac,gamma0*ones(size(omegac)),'wo','markersize',3);
    end
    daspect([2 1 1]);
    if bprint, fontsize(14); print(gcf,'-depsc',sprintf('%s/Vhat_time_all_%s',pdir,gtag)); end
    
    % histogram
    vtot = NaN(ngamma,1);
    for ii=1:ngamma
        vtemp = Vplot(ii,:);
        vtot(ii) = sum(vtemp(~isnan(vtemp)));
    end
    figure; plot(gammavec*deg,vtot/3600,'.-'); grid on;
    xlabel('\gamma, deg'); ylabel('time to compute Vhat, hours');
    if bprint, fontsize(14); print(gcf,'-depsc',sprintf('%s/Vhat_time_hist_%s',pdir,gtag)); end
end

%==========================================================================
