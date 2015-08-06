function plotMT_eigvec(pl1,az1,pl2,az2,pl3,az3)
%PLOTMT_EIGVEC plot three eigenvectors as trend-plunge angles
%
% See U2pa.m to convert from basis to plunge-azimuth (trend = azimuth).
%

% number of bins in each histogram
nbin = 18;      % 36
nedge = nbin + 1;
hflat = 1/nbin;

%plbin = linspace(0,90,nedge);
plbin = asin(linspace(0,1,nedge))*180/pi;
azbin = linspace(0,360,nedge);

ptick = [0:15:90];
atick = [0:60:360];

figure; nr=3; nc=2;
subplot(nr,nc,1); hold on; plot_histo(pl1,plbin); 
plot([0 90],hflat*[1 1],'r','linewidth',2);
xlabel('plunge angle of eigvec #1'); set(gca,'xtick',ptick);

subplot(nr,nc,2); hold on; plot_histo(az1,azbin);
plot([0 360],hflat*[1 1],'r','linewidth',2);
xlabel('azimuth angle of eigvec #1'); set(gca,'xtick',atick);

subplot(nr,nc,3); hold on; plot_histo(pl2,plbin);
plot([0 90],hflat*[1 1],'r','linewidth',2);
xlabel('plunge angle of eigvec #2'); set(gca,'xtick',ptick);

subplot(nr,nc,4); hold on; plot_histo(az2,azbin);
plot([0 360],hflat*[1 1],'r','linewidth',2);
xlabel('azimuth angle of eigvec #2'); set(gca,'xtick',atick);

subplot(nr,nc,5); hold on; plot_histo(pl3,plbin);
plot([0 90],hflat*[1 1],'r','linewidth',2);
xlabel('plunge angle of eigvec #3'); set(gca,'xtick',ptick);

subplot(nr,nc,6); hold on; plot_histo(az3,azbin);
plot([0 360],hflat*[1 1],'r','linewidth',2);
xlabel('azimuth angle of eigvec #3'); set(gca,'xtick',atick);

%==========================================================================
