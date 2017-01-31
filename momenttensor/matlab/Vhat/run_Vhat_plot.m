%
% run_Vhat_plot.m
%
% Compute a library of Vhat'(omega) curves that show the uniform
% distribution of moment tensors for a fixed source type (i.e., fixed lune
% point or fixed eigenvalues).
%
% W. Tape and C. Tape, GJI 2016
% A confidence parameter for seismic moment tensors.
%

close all
clc
clear

deg = 180/pi;
%==========================
bplot_all = false;
bprint = false;

path_Vhat;
bdir = Vhatdir_write;
pdir = Vhatdir_figs;
xticks = [0:30:180];
ax0 = [-1/deg 181/deg -0.02 1.1];
stvres = 'Vhat_{\gamma}(\omega)  -  Vhat_{-\gamma}(\omega)';
delta_deg = 0;
%==========================

% load V_gamma(omega) and V_gamma'(omega) for gamma > 0
USE_POSITIVE_GAMMA = true;
[pVp,pWp,pGp,pV,pW,pG] = load_Vgammaomega(USE_POSITIVE_GAMMA);
omegalib_rad = pW(1,:)';
omegalibp_rad = pWp(1,:)';
gamma_rad_pos = pG(:,1);

% cut gamma = 0 and gamma = 30, since the focus here is on the residual
% plots with negative gamma (gamma = -30 is not calculatable)
pVp([1 end],:) = [];
pWp([1 end],:) = [];
pGp([1 end],:) = [];
pV([1 end],:) = [];
pW([1 end],:) = [];
pG([1 end],:) = [];
gamma_rad_pos([1 end]) = [];

% load V_gamma(omega) and V_gamma'(omega) for gamma > 0
USE_POSITIVE_GAMMA = false;
[nVp,nWp,nGp,nV,nW,nG] = load_Vgammaomega(USE_POSITIVE_GAMMA);
gamma_rad_neg = nG(:,1);

% residual (note the convention)
Vres_all = pV - nV;

% pre-computed critical values
[omegacritcell,gamma_deg_crit] = load_omegacrit;
omegacritcell([1 end]) = [];

%----------------------------------------

omega = omegalib_rad;
gammavec = gamma_rad_pos;
nomega = length(omega);
ngamma = length(gammavec);

% CHANGE THESE INDICES AND USE bplot_all=true TO SEE INDIVIDUAL PLOTS
kmin = 1; kmax = ngamma;
%kmin = ngamma; kmax = kmin;
%kmin = 15; kmax = kmin;

for kk=kmin:kmax
    gamma_deg = gammavec(kk)*deg;
    omegacrit = omegacritcell{kk};
    omegacrit = omegacrit(:,1);
    
    pVhat  = pV(kk,:);
    pVhatp = pVp(kk,:);
    nVhat  = nV(kk,:);
    nVhatp = nVp(kk,:);
    
    % set all values omega > omega4 to NaN
    omega4 = omegacrit(end);
    inan   = find(omegalib_rad > omega4/deg);
    inanp  = find(omegalibp_rad > omega4/deg);
    %pVhatp(inanp) = NaN;
    %pVhat(inan) = NaN;
    %nVhatp(inanp) = NaN;
    %nVhat(inan) = NaN;
    Vres_all(kk,inan) = NaN;
    
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

                h1 = plot(omega,pVhat,'b-','linewidth',1);
                h2 = plot(omega,nVhat,'r--','linewidth',1);
                legend([h1 h2],'Vhat_{\gamma}(\omega)','Vhat_{-\gamma}(\omega)','location','northwest');
                title(sprintf('Vhat(\\omega) for lune point (\\gamma = %.1f^\\circ, \\delta = %.1f^\\circ)',gamma_deg,delta_deg));
                ylabel('Vhat(\omega)'); axis equal, axis(ax0);
                %ylabel('Vhat(\omega) and Vhat''(\omega)');
            else
                %ymx = 1e-4;
                [rmx,imx] = max(abs(Vres_all(kk,:)));
                if rmx > 0, vmx = rmx; else vmx = -rmx; end
                ymx = 1.1*rmx;
                [xmat, ymat] = vertlines(omegacrit/deg,-ymx,ymx);
                plot(xmat,ymat,'k--');
                ytxt = 0.9*ymx;
                for jj=1:4, text((omegacrit(jj)+dxx)/deg,ytxt,sprintf('\\omega_%i',jj)); end

                h1 = plot(omega,Vres_all(kk,:),'k.-');
                %h2 = plot(omega,tocvn,'r.');
                %h3 = plot(omega,tocvp,'c.');
                %legend([h1],'t(Vn)+t(Vp), sec','t(Vn), sec','t(Vp), sec','location','northwest');
                ylabel(stvres);
                title(sprintf('max difference is %.1e at \\omega = %.1f deg',rmx,omega(imx)*deg));
                plot(omega(imx),rmx,'ro');
                axis([ax0(1:2) [-1 1]*ymx]);
            end
            set(gca,'xtick',xticks/deg,'xticklabel',numvec2cell(xticks));
            xlabel('\omega, deg');
        end

        if bprint
            ftag = sprintf('Vhat_100gamma_%4.4i_diff',round(gamma_deg*100));
            orient tall; print(gcf,'-dpdf',sprintf('%s/%s',pdir,ftag));
            pause(1); close all;
        end
    end  % bplot
end

%%
Vplot = Vres_all; ymx = 5e-7;
figure; hold on;
pcolor(pW*deg,pG*deg,Vplot); shading flat; axis([0 180 0 30]);
xlabel('\omega, deg'); ylabel('\gamma, degrees');
caxis([-1 1]*ymx); colorbar
title(stvres); 
for ii=1:length(gammavec)
    gamma0 = abs(gammavec(ii))*deg;
    delta0 = 0;
    lam = lune2lam(gamma0,delta0);
    omegac = lam2omegacrit(lam);
    omegac = omegacritmod(omegac);
    omegac(1) = [];
    plot(omegac,gamma0*ones(size(omegac)),'ko','markersize',3);
end
daspect([2 1 1]);
if bprint
    print(gcf,'-depsc',sprintf('%s/Vhat_diff_all',pdir));
end

%break

%% remove NaN from the calculation
figure; plot(pW(isnan(Vplot))*deg,pG(isnan(Vplot))*deg,'.');
xlabel('\omega, deg'); ylabel('\gamma, deg'); daspect([2 1 1]);
vplot = Vplot(:);
vplot(isnan(vplot)) = [];
rms_all = rms(vplot);
% check
rms(vplot)
norm(vplot) / sqrt(length(vplot))
sqrt( sum(vplot.^2) / length(vplot) )

%% histogram of residuals
ymx = 1e-6; edges = linspace(-1e-6,1e-6,39);
figure; plot_histo(vplot,edges);
xlabel(stvres);
title(sprintf('RMS = %.2e [%i entries]',rms_all,length(vplot)));
if bprint
    fontsize(14); print(gcf,'-depsc',sprintf('%s/Vhat_diff_hist',pdir));
end

%% RMS calculation as a function of gamma
rms_gamma = NaN(ngamma,1);
for ii=1:length(gammavec)
    vhat = Vplot(ii,:);
    rms_gamma(ii) = rms(vhat(~isnan(vhat)));
end
figure; hold on; grid on;
plot(gammavec*deg,rms_gamma,'.-');
plot([0 30],[1 1]*rms_all,'r--');
xlabel('\gamma, deg'); ylabel(sprintf('RMS[ %s ]',stvres));
xlim([0 30]);
if bprint
    fontsize(14); print(gcf,'-depsc',sprintf('%s/Vhat_rms_gamma',pdir));
end

%% RMS calculation as a function of omega
rms_omega = NaN(ngamma,1);
for kk=1:length(omega)
    vhat = Vplot(:,kk);
    rms_omega(kk) = rms(vhat(~isnan(vhat)));
end
figure; hold on; grid on;
plot(omega*deg,rms_omega,'.-');
plot([0 180],[1 1]*rms_all,'r--');
xlabel('\omega, deg'); ylabel(sprintf('RMS[ %s ]',stvres));
xlim([-1 181]);
if bprint
    fontsize(14); print(gcf,'-depsc',sprintf('%s/Vhat_rms_omega',pdir));
end

%==========================================================================
